#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
    CREATE EXTENSION postgis_raster;
    CREATE SCHEMA crops;
    CREATE SCHEMA dssat;
    CREATE SCHEMA vic;
EOSQL
