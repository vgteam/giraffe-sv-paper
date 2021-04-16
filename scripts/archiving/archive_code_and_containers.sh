#!/usr/bin/env bash

# archive_code_and_containers.sh: Save source code commits and Docker containers referenced in this repository and the paper.

set -e
set -x

TEX_FILES=($HOME/build/giraffe-paper/main.tex $HOME/build/giraffe-paper/supplement.tex)
DEST_DIR=/nanopore/cgl/data/giraffe
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
    TARBALL_FILE="$(realpath "${TARBALL_DIR}/${TOOL_NAME}-${REF}.tar.gz")"
    
    if [[ ! -e "${TARBALL_FILE}" ]] ; then
        if [[ "${TOOL_NAME}" == "vg" && "${REF}" == v*.*.* ]] ; then
            # vg ships premade tarballs for real releases
            curl -sSL "https://github.com/vgteam/vg/releases/download/${REF}/vg-${REF}.tar.gz" > "${TARBALL_FILE}"
        else
            # Go make a tarball ourselves
            mkdir -p "${TOOL_DIR}"
            mkdir -p "${TARBALL_DIR}"
            rm -Rf "${CLONE_DIR}"
            git clone "${CLONE_URL}" "${CLONE_DIR}"
            (cd "${CLONE_DIR}" && git fetch --tags origin && git checkout "${REF}" && git submodule update --init --recursive)
            rm -Rf "${CLONE_DIR}/.git"
            find "${CLONE_DIR}" -name ".git" -exec rm -Rf "{}" \;
            # Compress with a nice relative path
            (cd "${CLONE_DIR}/.." && tar -czf "${TARBALL_FILE}" "$(basename "${CLONE_DIR}")")
            rm -Rf "${CLONE_DIR}"
        fi
    fi
    
    if [[ "${TOOL_NAME}" == "vg" && "${REF}" == v*.*.* && ! -e "${TARBALL_DIR}/vg" ]] ; then
        # vg will also ship static Linux x86_64 binaries for official releases.
        curl -sSL "https://github.com/vgteam/vg/releases/download/${REF}/vg" > "${TARBALL_DIR}/vg"
        chmod +x "${TARBALL_DIR}/vg"
    fi
}


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
    if [[ "${FILENAME}" == *archive_code_and_containers.sh ]] ; then
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
    echo "vg commit: ${VG_COMMIT}"
    archive_ref vg https://github.com/vgteam/vg.git "${VG_COMMIT}"
done
for TOIL_VG_COMMIT in $(printf "%s\n" "${TOIL_VG_COMMITS[@]}" | sort | uniq) ; do
    echo "toil-vg commit: ${TOIL_VG_COMMIT}"
    archive_ref toil-vg https://github.com/vgteam/toil-vg.git "${TOIL_VG_COMMIT}"
done
for TOIL_COMMIT in $(printf "%s\n" "${TOIL_COMMITS[@]}" | sort | uniq) ; do
    echo "toil commit: ${TOIL_COMMIT}"
    archive_ref toil https://github.com/DataBiosphere/toil.git "${TOIL_COMMIT}"
done

for VG_DOCKER in $(printf "%s\n" "${VG_DOCKERS[@]}" | sort | uniq) ; do
    echo "vg docker: ${VG_DOCKER}"
    archive_container vg "${VG_DOCKER}"
done
for TOIL_DOCKER in $(printf "%s\n" "${TOIL_DOCKERS[@]}" | sort | uniq) ; do
    echo "toil docker: ${TOIL_DOCKER}"
    archive_container toil "${TOIL_DOCKER}"
done

rm -Rf "${WORK_DIR}"


