# Giraffe Paper Data Archive

This archive contains the software developed for the manuscript *Genotyping common, large structural variations in 5,202 genomes using pangenomes, the Giraffe mapper, and the vg toolkit*. It is organized as follows:

* `code`: Contains source code
    * `giraffe-sv-paper`: Contains scripts used to manage the paper workflows and define the analyses, from [the repository here](https://github.com/vgteam/giraffe-sv-paper)
        * `giraffe-sv-paper.bundle`: a [Git bundle](https://git-scm.com/docs/git-bundle) containing the history of the scripts up to the latest version when the archive was created
    * `toil`: Contains the [Toil workflow execution system](https://github.com/DataBiosphere/toil)
        * `<commit>` or `<version>`: Contains a version of Toil used in the manuscript
            * `toil-<commit or version>.tar.gz`: contains an archive of the contents of the Toil repository at that version
        * `toil.bundle`: a [Git bundle](https://git-scm.com/docs/git-bundle) containing the history of the Toil codebase up to the latest version when the archive was created
    * `toil-vg`: Contains the [toil-vg workflow](https://github.com/vgteam/toil-vg)
        * `<commit>`: Contains a version of toil-vg used in the manuscript
            * `toil-vg-<commit>.tar.gz`: contains an archive of the contents of the Toil repository at that version
        * `toil-vg.bundle`: a [Git bundle](https://git-scm.com/docs/git-bundle) containing the history of the toil-vg codebase up to the latest version when the archive was created
    * `vg`: Contains the [vg tool](https://github.com/vgteam/vg), which includes the Giraffe mapper
        * `<commit>`: Contains a version of vg used in the manuscript
            * `vg-<commit>.tar.gz`: contains an archive of the contents of the vg repository and all submodules at that commit
        * `<version>`: Contains a released version of vg used in the manuscript
            * `vg-<version>.tar.gz`: contains an archive of the contents of the vg repository and all submodules at that version
            * `vg`: A statically linked x86_64 Linux executable of this release of vg
        * `vg.bundle`: a [Git bundle](https://git-scm.com/docs/git-bundle) containing the history of the vg codebase up to the latest version when the archive was created. Does not include the contents of submodules.
* `containers`: Contains Docker containers, organized by project and then Docker tag. Source repository information is not included as tags are unique across the manuscript's workflows.
    * `toil`: Docker containers of Toil versions used in the manuscript
        * `<tag>.tar.gz`: Container of version identified by `<tag>`, exported via docker export` and compressed.
    * `vg`: Docker containers of vg used in the manuscript.
        * `<tag>.tar.gz`: Container of version identified by `<tag>`, exported via docker export` and compressed. Note that tag `giraffe-paper` corresponds to vg v1.31.0, commit 08faee067037ece539a237a008bcdefc84b681b0.


