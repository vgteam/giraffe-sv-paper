#!/usr/bin/env bash

# save_products.sh: Save work products that ought to be re-usable

set -x

DEST_DIR=/nanopore/cgl/data/giraffe

function web_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function download_if_exists() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}" || true
    fi
}

PRODUCTS_DIR="${DEST_DIR}/products"
mkdir -p "${PRODUCTS_DIR}"



# TODO: add lines here to actually save the product files:
#web_download https://google-storage.gov/whatever.dat "${PRODUCTS_DIR}/whatever.dat"
# Also document each in products-readme.md

download s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.full.gbwt "${PRODUCTS_DIR}/HGSVC_hs38d1.full.gbwt"
download s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.dist "${PRODUCTS_DIR}/HGSVC_hs38d1.dist"
download s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.xg "${PRODUCTS_DIR}/HGSVC_hs38d1.xg"

## SV-related files
SVFILES="vggiraffe-sv-eqtl-geuvadis.FDR01.csv
vggiraffe-geuvadis-sveqtl-gene-families.csv
vggiraffe-geuvadis-eqtl-svonly.csv
vggiraffe-geuvadis-eqtl-snv-indel-svs.csv.gz
vggiraffe-sv-superpop-af-diff-med10.csv.gz
vggiraffe-sv-2504kgp-pcgenes.tsv.gz
vggiraffe-sv-mesa-svsites.vcf.gz
vggiraffe-sv-mesa-svsites.vcf.gz.tbi
vggiraffe-sv-2504kgp-svsites.vcf.gz
vggiraffe-sv-2504kgp-svsites.vcf.gz.tbi
vggiraffe-sv-2504kgp-svsites.gt.vcf.gz
vggiraffe-sv-2504kgp-svsites.gt.vcf.gz.tbi
vggiraffe-sv-2504kgp-raw.vcf.gz
vggiraffe-sv-2504kgp-raw.vcf.gz.tbi
vggiraffe-sv-relkgp-raw.vcf.gz
vggiraffe-sv-relkgp-raw.vcf.gz.tbi"
for ff in $SVFILES ; do
    download s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/$ff "${PRODUCTS_DIR}/${ff}"
done


# Put the README in place
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
# Fix the links to the archive data to be relative
cat "${SCRIPT_DIR}/products-readme.md" | sed 's_https://cgl.gi.ucsc.edu/data/giraffe/__g' > "${PRODUCTS_DIR}/README.md"

for PRODUCT_FILENAME in $(ls "${PRODUCTS_DIR}") ; do
    # Upload each product to Zenodo
    if [[ ! -z "${ZENODO_DEPOSITION}" && ! -z "${ZENODO_TOKEN}" ]] ; then
        export FILEPATH="${PRODUCTS_DIR}/${PRODUCT_FILENAME}"
        # Upload the product file onto the Zenodo deposition specified by the environment
        python3 -c 'import requests; 
import os;
deposition=os.environ["ZENODO_DEPOSITION"]; 
filepath=os.environ["FILEPATH"]; 
filename=os.path.basename(filepath); 
params={"access_token": os.environ["ZENODO_TOKEN"]}; 
bucket=requests.get(f"https://www.zenodo.org/api/deposit/depositions/{deposition}", params=params).json()["links"]["bucket"]; 
requests.put(f"{bucket}/products/{filename}", data=open(filepath, "rb"), params=params).raise_for_status();'
    fi
done

chmod -R g+rw "${DEST_DIR}" 2>/dev/null || true
