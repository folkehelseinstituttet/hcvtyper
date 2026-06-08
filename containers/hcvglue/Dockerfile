# hcvglue-allinone: single hermetic container for per-BAM HCV-GLUE genotyping.
#
# Collapses the legacy 2-container `cvrbioinformatics/gluetools` +
# `cvrbioinformatics/gluetools-mysql` `--link` setup into one image:
# the GLUE engine (Java + BLAST + MAFFT) already lives in the base, and we add a
# MySQL 5.7 server + a patched config pointing the engine at a task-local mysqld
# (127.0.0.1). The HCV project SQL dump is NOT baked in (D-03) — it is staged at
# runtime via params.hcvglue_db and loaded by run-glue.sh.
#
# Publish target: ghcr.io/folkehelseinstituttet/hcvglue-allinone
# Build:  docker build -t ghcr.io/folkehelseinstituttet/hcvglue-allinone:1.1.114 containers/hcvglue/
#
# SECURITY: pin the base image by digest once verified, e.g.
#   FROM cvrbioinformatics/gluetools:1.1.114@sha256:<digest>
# and tag+digest-pin this published image in modules/local/hcvglue/main.nf
# (RESEARCH Open Question 1 / threat T-05-02).
FROM cvrbioinformatics/gluetools:1.1.114

# Add MySQL 5.7 server at IMAGE BUILD TIME (on the maintainer's build host),
# NOT inside a pipeline container at runtime — this does not violate GLUE-06.
#
# FALLBACK (RESEARCH A1): if the engine base lacks a usable apt source for
# mysql-server-5.7, replace the RUN below with a multi-stage copy of the mysqld
# binaries from the official mysql:5.7.36 image, e.g.:
#   COPY --from=mysql:5.7.36 /usr/sbin/mysqld /usr/sbin/mysqld
#   COPY --from=mysql:5.7.36 /usr/bin/mysql /usr/bin/mysqladmin /usr/bin/
# Documented here as the sanctioned fallback path.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends mysql-server-5.7 \
    && rm -rf /var/lib/apt/lists/*

# Patch the GLUE engine DB config: jdbcUrl host gluetools-mysql -> 127.0.0.1.
COPY gluetools-config.xml ${GLUE_HOME}/conf/gluetools-config.xml

# Bake the hermetic per-BAM entrypoint (MySQL lifecycle + GLUE report commands).
COPY run-glue.sh /usr/local/bin/run-glue.sh
RUN chmod +x /usr/local/bin/run-glue.sh
