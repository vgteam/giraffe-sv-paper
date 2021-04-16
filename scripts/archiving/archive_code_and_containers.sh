#!/usr/bin/env bash

# archive_code_and_containers.sh: Save source code commits and Docker containers referenced in this repository and the paper.

set -e
#set -x

TEX_FILES=($HOME/build/giraffe-paper/main.tex $HOME/build/giraffe-paper/supplement.tex)
DEST_DIR=/nanopore/cgl/data/giraffe/code


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
done
for TOIL_VG_COMMIT in $(printf "%s\n" "${TOIL_VG_COMMITS[@]}" | sort | uniq) ; do
    echo "toil-vg commit: ${TOIL_VG_COMMIT}"
done
for TOIL_COMMIT in $(printf "%s\n" "${TOIL_COMMITS[@]}" | sort | uniq) ; do
    echo "toil commit: ${TOIL_COMMIT}"
done

for VG_DOCKER in $(printf "%s\n" "${VG_DOCKERS[@]}" | sort | uniq) ; do
    echo "vg docker: ${VG_DOCKER}"
done
for TOIL_DOCKER in $(printf "%s\n" "${TOIL_DOCKERS[@]}" | sort | uniq) ; do
    echo "toil docker: ${TOIL_DOCKER}"
done


