#!/usr/bin/env bash

# archive_code_and_containers.sh: Save source code commits and Docker containers referenced in this repository and the paper.

set -e
set -x

TEX_FILES=($HOME/build/giraffe-paper/main.tex $HOME/build/giraffe-paper/supplement.tex)
BASE_DEST_DIR=/nanopore/cgl/data/giraffe
DEST_DIR=${BASE_DEST_DIR}/software
WORK_DIR="$(mktemp -d)"

function archive_container {
    # Save a Docker container
    # archive_container TOOL_NAME CONTAINER_SPEC
    TOOL_NAME="${1}"
    CONTAINER_SPEC="${2}"
    TAG="$(echo "${CONTAINER_SPEC}" | cut -f2 -d':')"
    
    CONTAINER_DIR="${DEST_DIR}/containers/${TOOL_NAME}"
    CONTAINER_TAR="${WORK_DIR}/${TAG}.tar"
    CONTAINER_FILE="${CONTAINER_DIR}/${TAG}.tar.gz"
    
    if [[ ! -e "${CONTAINER_FILE}" ]] ; then
        mkdir -p "${CONTAINER_DIR}"
        docker pull "${CONTAINER_SPEC}"
        docker save "${CONTAINER_SPEC}" -o "${CONTAINER_TAR}"
        pigz "${CONTAINER_TAR}"
        mv "${CONTAINER_TAR}.gz" "${CONTAINER_FILE}"
    fi
    chmod 644 "${CONTAINER_FILE}"
}

function archive_ref {
    # Save a ref of a Git repository
    # archive_ref TOOL_NAME CLONE_URL REF
    TOOL_NAME="${1}"
    CLONE_URL="${2}"
    REF="${3}"
    
    TOOL_DIR="${DEST_DIR}/code/${TOOL_NAME}"
    CLONE_DIR="${WORK_DIR}/${TOOL_NAME}-${REF}"
    TARBALL_DIR="${TOOL_DIR}/${REF}"
    TARBALL_FILE="${TARBALL_DIR}/${TOOL_NAME}-${REF}.tar.gz"
    
    if [[ ! -e "${TARBALL_FILE}" ]] ; then
        mkdir -p "${TARBALL_DIR}"
        if [[ "${TOOL_NAME}" == "vg" && "${REF}" == v*.*.* ]] ; then
            # vg ships premade tarballs for real releases
            curl -sSL "https://github.com/vgteam/vg/releases/download/${REF}/vg-${REF}.tar.gz" > "${TARBALL_FILE}"
        else
            # Go make a tarball ourselves
            rm -Rf "${CLONE_DIR}"
            git clone "${CLONE_URL}" "${CLONE_DIR}"
            (cd "${CLONE_DIR}" && git fetch --tags origin && (git checkout "${REF}" || git checkout "releases/${REF}") && git submodule update --init --recursive)
            rm -Rf "${CLONE_DIR}/.git"
            find "${CLONE_DIR}" -name ".git" -exec rm -Rf "{}" \;
            # Compress with a nice relative path
            TARBALL_ABSPATH="$(realpath "${TARBALL_FILE}")"
            (cd "${CLONE_DIR}/.." && tar -czf "${TARBALL_ABSPATH}" "$(basename "${CLONE_DIR}")")
            rm -Rf "${CLONE_DIR}"
        fi
    fi
    chmod 644 "${TARBALL_FILE}"
    
    if [[ "${TOOL_NAME}" == "vg" && "${REF}" == v*.*.* ]] ; then
        # vg will also ship static Linux x86_64 binaries for official releases.
        if [[ ! -e "${TARBALL_DIR}/vg" ]] ; then
            curl -sSL "https://github.com/vgteam/vg/releases/download/${REF}/vg" > "${TARBALL_DIR}/vg"
        fi
        chmod 755 "${TARBALL_DIR}/vg"
    fi
}

function bundle_all {
    # Save all refs of a Git repository as a bundle
    # bundle_refs TOOL_NAME CLONE_URL
    # Can't really work with submodules.
    TOOL_NAME="${1}"
    shift
    CLONE_URL="${1}"
    shift
    
    TOOL_DIR="${DEST_DIR}/code/${TOOL_NAME}"
    CLONE_DIR="${WORK_DIR}/${TOOL_NAME}"
    BUNDLE_DIR="${TOOL_DIR}"
    BUNDLE_FILE="${BUNDLE_DIR}/${TOOL_NAME}.bundle"
    
    mkdir -p "${BUNDLE_DIR}"
    rm -Rf "${CLONE_DIR}"
    git clone --mirror "${CLONE_URL}" "${CLONE_DIR}"
    (cd "${CLONE_DIR}" && git fetch --tags origin)
    
    BUNDLE_ABSPATH="$(realpath "${BUNDLE_FILE}")"
    (cd "${CLONE_DIR}" && git bundle create "${BUNDLE_ABSPATH}" --all)
    rm -Rf "${CLONE_DIR}"
    
    chmod 644 "${BUNDLE_FILE}"
}

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

PRODUCTS_DIR="${BASE_DEST_DIR}/products"

SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

ALL_DOCKERS=()
ALL_GIT_URLS=()

VG_DOCKERS=()
VG_COMMITS=()
TOIL_VG_GIT_URLS=()
TOIL_VG_COMMITS=()
TOIL_DOCKERS=()
TOIL_COMMITS=()

for FILENAME in "${TEX_FILES[@]}" ; do
    for DOCKER in $(cat "${FILENAME}" | grep -v "^ *%" | grep -o -E '\\docker{[^}]*}' | sed -e 's/^.*{//g' -e 's/}$//g') ; do
        ALL_DOCKERS+=("${DOCKER}")
    done
    for COMMIT in $(cat "${FILENAME}" | grep -v "^ *%" | grep -o -E '\\vgcommit{[^}]*}' | sed -e 's/^.*{//g' -e 's/}$//g') ; do
        VG_COMMITS+=("${COMMIT}")
    done
    for COMMIT in $(cat "${FILENAME}" | grep -v "^ *%" | grep -o -E '\\toilvgcommit{[^}]*}' | sed -e 's/^.*{//g' -e 's/}$//g') ; do
        TOIL_VG_COMMITS+=("${COMMIT}")
    done
done

for FILENAME in $(find "${SCRIPT_DIR}/..") ; do
    if [[ -d "${FILENAME}" ]] ; then
        # Skip directories
        continue
    fi
    if [[ "${FILENAME}" == *archive_code_containers_and_products.sh ]] ; then
        # Skip ourselves
        continue
    fi
    for DOCKER in $(cat "${FILENAME}" | grep -o -E "[0-9a-zA-Z./_-]+/(vg|toil)(:[0-9a-zA-Z._-]+)") ; do
        ALL_DOCKERS+=("${DOCKER}")
    done
    for URL in $(cat "${FILENAME}" | grep -o -E '[0-9a-zA-Z+_/.:-]+\.git@[0-9a-zA-Z+_#=-]+') ; do
        ALL_GIT_URLS+=("${URL}")
    done
    for VERSION in $(cat "${FILENAME}" | grep -o -E 'vg v[0-9]+\.[0-9]+\.[0-9]+' | cut -f2 -d' ') ; do
        VG_COMMITS+=("${VERSION}")
    done
done

for DOCKER in "${ALL_DOCKERS[@]}" ; do
    if [[ "${DOCKER}" == quay.io/vgteam/vg:* ]] ; then
        VG_DOCKERS+=("${DOCKER}")
        if [[ "${DOCKER}" == quay.io/vgteam/vg:ci-* ]] ; then
            VG_COMMITS+=("$(echo "${DOCKER}" | cut -f2 -d':' | cut -f3 -d'-')")
        elif [[ "${DOCKER}" == quay.io/vgteam/vg:v*.*.* ]] ; then
            # Use the version tag as a commit
            VG_COMMITS+=("$(echo "${DOCKER}" | cut -f2 -d':')")
        fi
    elif [[ "${DOCKER}" == xhchang/vg:* ]] ; then
        VG_DOCKERS+=("${DOCKER}")
    elif [[ "${DOCKER}" == quay.io/ucsc_cgl/toil:* ]] ; then
        TOIL_DOCKERS+=("${DOCKER}")
        if [[ "${DOCKER}" == quay.io/ucsc_cgl/toil:*-*-* ]] ; then
            # Toil dev version
            TOIL_COMMITS+=("$(echo "${DOCKER}" | cut -f2 -d':' | cut -f2 -d'-')")
        elif [[ "${DOCKER}" == quay.io/ucsc_cgl/toil:*-* ]] ; then
            # Toil release
            TOIL_COMMITS+=("$(echo "${DOCKER}" | cut -f2 -d':' | cut -f1 -d'-')")
        fi
    else
        echo "Unrecognized Docker: ${DOCKER}"
    fi
done

for URL in "${ALL_GIT_URLS[@]}" ; do
    if [[ "${URL}" == git+https://github.com/vgteam/toil-vg.git@* ]] ; then
        TOIL_VG_GIT_URLS+=("${URL}")
        TOIL_VG_COMMITS+=("$(echo "${URL}" | cut -f2 -d'@' | cut -f1 -d'#')")
    else
        echo "Unrecognized Git URL: ${URL}"
    fi
done

for VG_COMMIT in $(printf "%s\n" "${VG_COMMITS[@]}" | sort | uniq) ; do
    if [[ "${VG_COMMIT}" == "file" ]] ; then
        # Not sure how this go in
        continue
    fi
    echo "vg commit: ${VG_COMMIT}"
    archive_ref vg https://github.com/vgteam/vg.git "${VG_COMMIT}"
done
bundle_all vg https://github.com/vgteam/vg.git

for TOIL_VG_COMMIT in $(printf "%s\n" "${TOIL_VG_COMMITS[@]}" | sort | uniq) ; do
    echo "toil-vg commit: ${TOIL_VG_COMMIT}"
    archive_ref toil-vg https://github.com/vgteam/toil-vg.git "${TOIL_VG_COMMIT}"
done
bundle_all toil-vg https://github.com/vgteam/toil-vg.git

for TOIL_COMMIT in $(printf "%s\n" "${TOIL_COMMITS[@]}" | sort | uniq) ; do
    echo "toil commit: ${TOIL_COMMIT}"
    archive_ref toil https://github.com/DataBiosphere/toil.git "${TOIL_COMMIT}"
done
bundle_all toil https://github.com/DataBiosphere/toil.git

for VG_DOCKER in $(printf "%s\n" "${VG_DOCKERS[@]}" | sort | uniq) ; do
    echo "vg docker: ${VG_DOCKER}"
    archive_container vg "${VG_DOCKER}"
done
for TOIL_DOCKER in $(printf "%s\n" "${TOIL_DOCKERS[@]}" | sort | uniq) ; do
    echo "toil docker: ${TOIL_DOCKER}"
    archive_container toil "${TOIL_DOCKER}"
done

# Now archive all the paper scripts and history.
bundle_all giraffe-sv-paper ${SCRIPT_DIR}/../..

# And put the software README in place
cp "${SCRIPT_DIR}/software-readme.md" "${DEST_DIR}/README.md"

rm -Rf "${WORK_DIR}"

chmod -R g+rw "${DEST_DIR}" 2>/dev/null || true

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
# Fix the links to the archive data to be relative
cat "${SCRIPT_DIR}/products-readme.md" | sed 's_https://cgl.gi.ucsc.edu/data/giraffe/__g' > "${PRODUCTS_DIR}/README.md"

# Zip all the software and products together. Use zip because incremental update of touched
# files is efficient. Do it all at once because Zenodo seems to refuse the second zip otherwise.
COMBINED_ZIP_FILE="archive.zip"
COMBINED_ZIP_ABSPATH="$(realpath "${COMBINED_ZIP_FILE}")"
(cd "${BASE_DEST_DIR}" && zip -ur "${COMBINED_ZIP_ABSPATH}" "$(basename "${DEST_DIR}")" "$(basename "${PRODUCTS_DIR}")")

chmod  g+rw "${COMBINED_ZIP_FILE}" 2>/dev/null || true

if [[ ! -z "${ZENODO_DEPOSITION}" && ! -z "${ZENODO_TOKEN}" ]] ; then
    export FILEPATH="${COMBINED_ZIP_FILE}"
    python3 -c 'import requests; 
import os;
deposition=os.environ["ZENODO_DEPOSITION"]; 
filepath=os.environ["FILEPATH"]; 
filename=os.path.basename(filepath); 
params={"access_token": os.environ["ZENODO_TOKEN"]}; 
bucket=requests.get(f"https://www.zenodo.org/api/deposit/depositions/{deposition}", params=params).json()["links"]["bucket"]; 
requests.put(f"{bucket}/{filename}", data=open(filepath, "rb"), params=params).raise_for_status();'
fi

