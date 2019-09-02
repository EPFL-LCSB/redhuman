redHUMAN
========

Reduction of human genome-scale models. Paper: Maria Masid, Meriç Ataman and Vassily Hatzimanikatis. "redHUMAN: analyzing human metabolism and growth media through systematic reductions of thermodynamically curated genome-scale models" 

Requirements
------------

You will need to have `Git LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/EPFL-LCSB/redhuman.git /path/to/redhuman
    cd /path/to/redhuman
    git lfs install
    git lfs pull

The scripts have been developed with Matlab 2017b, and CPLEX 12.7 (freely downloadable with the `IBM Academic initiative <https://developer.ibm.com/academic/>`_), and successfully ran on several other versions of both softwares. However, it is important to respect the IBM compatibility specs sheets between Matlab, CPLEX, and the computer OS - available `on IBM's website <https://www.ibm.com/software/reports/compatibility/clarity/index.html>`_.

This module requires `matTFA <https://github.com/EPFL-LCSB/mattfa/>`_ and `redGEM <https://github.com/EPFL-LCSB/redgem/>`_.

Generating reduced models
-------------------------
1. Place the corresponding thermodynamic data from the `redhuman data <https://github.com/EPFL-LCSB/redhuman/data>`_ folder into the `matTFA thermoDatabases <https://github.com/EPFL-LCSB/matTFA/thermoDatabases>`_ folder.
2. Place the corresponding curated GEM from the `redhuman GEMs <https://github.com/EPFL-LCSB/redhuman/redhuman/GEMs>`_ folder into the `redGEM GEMs <https://github.com/EPFL-LCSB/redgem/GEMs>`_ folder.
3. Place the corresponding get_redHUMAN file from the `redhuman <https://github.com/EPFL-LCSB/redhuman/redhuman>`_ folder into the `redGEM runFileExample <https://github.com/EPFL-LCSB/redgem/runFiles>`_  folder.
4. Modify the get_redHUMAN file to add the paths and the dessired parameters
5. Run the get_redHUMAN file

Model validation
----------------
Run the scripts from the `postprocessing <https://github.com/EPFL-LCSB/redhuman/redhuman/postprocessing>`_ folder to test the metabolic tasks, the gene essentiality analysis and the flux variability analysis as they are done for the paper.


License
=======
The software in this repository is put under an APACHE licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/redhuman/blob/master/LICENSE>`_ file for more details.

