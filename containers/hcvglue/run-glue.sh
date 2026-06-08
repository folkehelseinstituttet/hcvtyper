#!/bin/bash
#
# run-glue.sh - Hermetic per-BAM HCV-GLUE entrypoint for the hcvglue-allinone image.
#
# Replaces the legacy bin/run_hcvglue.sh 2-container linked-MySQL dance. This script
# runs entirely INSIDE a single container task: it initialises a task-local
# MySQL 5.7 datadir, starts mysqld bound to loopback, loads the staged HCV-GLUE
# SQL dump, runs the GLUE reportBam / reportBamAsHtml commands, then shuts mysqld
# down cleanly. It spawns no sibling containers and touches no host Docker socket.
#
# Positional args (matches the module script invocation):
#   $1 = DUMP       staged hcvglue_db (.sql.gz), e.g. ncbi_hcv_glue.sql.gz
#   $2 = THRESHOLD  params.hcvglue_threshold
#   $3 = BAM        the single staged BAM (basename drives output names)
#
# Sources:
#   - github.com/giffordlabcvr/gluetools installGlueProject.sh (gunzip|mysql load)
#   - github.com/giffordlabcvr/gluetools init_glue_db.sql (GLUE_TOOLS + gluetools user)
#   - bin/run_hcvglue.sh:143-169 (preserved reportBam / reportBamAsHtml command shape)
#   - Singularity-safe task-writable datadir override (Apptainer read-only rootfs)
#
set -euo pipefail

DUMP="$1"       # staged hcvglue_db (.sql.gz)
THRESHOLD="$2"  # params.hcvglue_threshold
BAM="$3"        # single staged BAM

# 1. Task-writable MySQL locations (CRITICAL for Singularity read-only rootfs +
#    non-root invoking UID). NEVER use the image default system datadir.
DATADIR="$PWD/mysql-data"
SOCK="$PWD/mysqld.sock"
PIDF="$PWD/mysqld.pid"
MYSQL_USER="$(id -un)"

mkdir -p "$DATADIR"

# Initialise a fresh datadir. mysql 5.7 uses --initialize-insecure; fall back to
# the older mysql_install_db for builds that ship it instead.
mysqld --initialize-insecure --datadir="$DATADIR" --user="$MYSQL_USER" 2>/dev/null \
  || mysql_install_db --datadir="$DATADIR" --user="$MYSQL_USER"

# 2. Start mysqld bound to loopback TCP 3306 (GLUE connects via JDBC TCP per the
#    patched gluetools-config.xml -> 127.0.0.1). Do NOT use --skip-networking.
mysqld --datadir="$DATADIR" --socket="$SOCK" --pid-file="$PIDF" \
       --bind-address=127.0.0.1 --port=3306 &
MYSQLD_PID=$!

# 3. Deterministic readiness wait (not a sleep guess), up to 120s.
ready=0
for _ in $(seq 1 120); do
  if mysqladmin --socket="$SOCK" ping 2>/dev/null | grep -q "is alive"; then
    ready=1
    break
  fi
  sleep 1
done
if [ "$ready" -ne 1 ]; then
  echo "ERROR: mysqld did not become ready within 120s" >&2
  kill "$MYSQLD_PID" 2>/dev/null || true
  exit 1
fi

# 4. Bootstrap GLUE_TOOLS db + gluetools user (from upstream init_glue_db.sql),
#    creating the user for both 127.0.0.1 and localhost so the 5.7 auth plugin
#    accepts the loopback JDBC connection, then load the dump.
mysql --socket="$SOCK" -u root <<'SQL'
CREATE DATABASE IF NOT EXISTS GLUE_TOOLS CHARACTER SET UTF8;
CREATE USER IF NOT EXISTS 'gluetools'@'127.0.0.1' IDENTIFIED BY 'glue12345';
CREATE USER IF NOT EXISTS 'gluetools'@'localhost' IDENTIFIED BY 'glue12345';
GRANT ALL ON GLUE_TOOLS.* TO 'gluetools'@'127.0.0.1';
GRANT ALL ON GLUE_TOOLS.* TO 'gluetools'@'localhost';
FLUSH PRIVILEGES;
SQL

gunzip -c "$DUMP" | mysql --socket="$SOCK" -u gluetools -pglue12345 GLUE_TOOLS

# 5. Run GLUE. The engine reads the patched gluetools-config.xml (jdbcUrl host =
#    127.0.0.1). Command shape preserved from bin/run_hcvglue.sh:143-169.
#    `|| true` keeps a single bad BAM from failing the whole task (matches legacy).
gluetools.sh -p cmd-result-format:json -EC \
  -i project hcv module phdrReportingController invoke-function reportBam "$BAM" "$THRESHOLD" \
  > "${BAM%.bam}.json" || true

gluetools.sh --console-option log-level:FINEST \
  --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml "$BAM" "$THRESHOLD" "${BAM%.bam}.html" || true

# 6. Clean shutdown of the task-local mysqld (kill fallback if mysqladmin fails).
mysqladmin --socket="$SOCK" -u root shutdown 2>/dev/null || kill "$MYSQLD_PID" 2>/dev/null || true
wait "$MYSQLD_PID" 2>/dev/null || true

echo "HCV-GLUE per-BAM analysis completed for $BAM"
