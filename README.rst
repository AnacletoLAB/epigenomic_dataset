epigenomic_dataset
=========================================================================================
|travis| |sonar_quality| |sonar_maintainability|
|codacy| |code_climate_maintainability| |pip| |downloads|

Python package wrapping ENCODE epigenomic data
for several reference cell lines.

How do I install this package?
----------------------------------------------
As usual, just download it using pip:

.. code:: shell

    pip install epigenomic_dataset

Tests Coverage
----------------------------------------------
Since some software handling coverages sometimes get slightly
different results, here's three of them:

|coveralls| |sonar_coverage| |code_climate_coverage|


TODO: THE FOLLOWING SECTION WILL NEED RESTRUCTURING IN A LITTLE BIT!

Preprocessed data for cis-regulatory regions
-----------------------------------------------
We have already downloaded and obtained the max window value for each promoter and enhancer
region for the cell lines A549, GM12878, H1, HEK293, HepG2, K562 and MCF7 in the dataset Fantom
and cell lines A549, GM12878, H1, HepG2 and K562 for the Roadmap dataset taking in consideration
all the target features listed `in the complete table of epigenomes <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/epigenomic_dataset/epigenomes.csv>`__.

The thresholds used for classifying the activations of enhancers and promoters in Fantom are the
default explained in the sister pipeline `CRR labels <https://github.com/LucaCappelletti94/crr_labels>`__
which handles the download and preprocessing of the data from Fantom and Roadmap.

=========  ==========  =============  =========  ===========  ==================================================================================================================================================
Dataset    Assembly      Window Size  Region     Cell line    Download URL
=========  ==========  =============  =========  ===========  ==================================================================================================================================================
fantom     hg38                  256  promoters  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/GM12878.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/A549.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/HEK293.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/HepG2.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/K562.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/H1.csv.xz?raw=true>`__
fantom     hg38                  256  promoters  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/promoters/MCF-7.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/GM12878.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/A549.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/HEK293.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/HepG2.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/K562.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/H1.csv.xz?raw=true>`__
fantom     hg38                  256  enhancers  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/256/enhancers/MCF-7.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/GM12878.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/A549.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/HEK293.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/HepG2.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/K562.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/H1.csv.xz?raw=true>`__
fantom     hg38                  128  promoters  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/promoters/MCF-7.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/GM12878.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/A549.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/HEK293.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/HepG2.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/K562.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/H1.csv.xz?raw=true>`__
fantom     hg38                  128  enhancers  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/128/enhancers/MCF-7.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/GM12878.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/A549.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/HEK293.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/HepG2.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/K562.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/H1.csv.xz?raw=true>`__
fantom     hg38                   64  promoters  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/promoters/MCF-7.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/GM12878.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/A549.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/HEK293.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/HepG2.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/K562.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/H1.csv.xz?raw=true>`__
fantom     hg38                   64  enhancers  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/64/enhancers/MCF-7.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/GM12878.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/A549.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/HEK293.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/HepG2.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/K562.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/H1.csv.xz?raw=true>`__
fantom     hg38                 1024  promoters  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/promoters/MCF-7.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/GM12878.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/A549.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/HEK293.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/HepG2.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/K562.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/H1.csv.xz?raw=true>`__
fantom     hg38                 1024  enhancers  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/1024/enhancers/MCF-7.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/GM12878.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/A549.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/HEK293.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/HepG2.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/K562.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/H1.csv.xz?raw=true>`__
fantom     hg38                  512  promoters  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/promoters/MCF-7.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  GM12878      `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/GM12878.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  A549         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/A549.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  HEK293       `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/HEK293.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  HepG2        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/HepG2.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  K562         `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/K562.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  H1           `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/H1.csv.xz?raw=true>`__
fantom     hg38                  512  enhancers  MCF-7        `Download <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/hg38/512/enhancers/MCF-7.csv.xz?raw=true>`__
=========  ==========  =============  =========  ===========  ==================================================================================================================================================

Here are the labels for all the considered cell lines.

+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+
|   Dataset         |   Promoters                                                                                                                                                                                                                                                           |   Enhancers                                                                                                                                                                                                                                                           |
+===================+==================================================================================================================================+====================================================================================================================================+==================================================================================================================================+====================================================================================================================================+
| Fantom            | `200 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/200/promoters.bed.gz?raw=true>`__  | `1000 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/1000/promoters.bed.gz?raw=true>`__  | `200 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/200/enhancers.bed.gz?raw=true>`__  | `1000 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/fantom/1000/enhancers.bed.gz?raw=true>`__  |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+
| Roadmap           | `200 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/roadmap/200/promoters.bed.gz?raw=true>`__ | `1000 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/roadmap/1000/promoters.bed.gz?raw=true>`__ | `200 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/roadmap/200/enhancers.bed.gz?raw=true>`__ | `1000 <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed/roadmap/1000/enhancers.bed.gz?raw=true>`__ |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+

TODO: align promoters and enhancers in a reference labels dataset.

The complete pipeline used to retrieve the CRR epigenomic data is available
`here <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/run_crr_build.py>`__.

Automatic retrieval of preprocessed data
----------------------------------------------
You can automatically retrieve the data as follows:

.. code:: python

    from epigenomic_dataset import load_epigenomes

    X, y = load_epigenomes(
        cell_line = "K562",
        dataset = "fantom",
        regions = "promoters",
        window_size = 200,
        root = "datasets" # Path where to download data
    )

Pipeline for epigenomic data
----------------------------------------------
The considered raw data are from `this query from the ENCODE project <https://www.encodeproject.org/search/?searchTerm=fold+change+over+control&type=Experiment&assembly=hg19&status=released&biosample_ontology.classification=cell+line&files.file_type=bigWig&replication_type=isogenic&audit.ERROR.category%21=extremely+low+read+depth&audit.ERROR.category%21=inconsistent+genetic+modification+reagent+source+and+identifier&audit.ERROR.category%21=missing+control+alignments&audit.ERROR.category%21=extremely+low+read+length&audit.NOT_COMPLIANT.category%21=insufficient+read+depth&audit.NOT_COMPLIANT.category%21=missing+controlled_by&audit.NOT_COMPLIANT.category%21=insufficient+read+length&audit.NOT_COMPLIANT.category%21=insufficient+replicate+concordance&audit.NOT_COMPLIANT.category%21=severe+bottlenecking&audit.NOT_COMPLIANT.category%21=control+insufficient+read+depth&audit.NOT_COMPLIANT.category%21=poor+library+complexity&limit=all>`_

You can find the `complete table of the available epigenomes here <https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/epigenomic_dataset/epigenomes.csv>`_.
These datasets were selected to have
(at time of the writing,  07/02/2020)
the least possible amount of known problems, such as
low read resolution.

You can run the pipeline as follows: suppose you
want to extract the epigenomic features for the cell lines HepG2 and H1:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"]
    )

If you want to specify where to store the files use:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"],
        path="path/to/my/target"
    )

By default, the downloaded bigWig files are not deleted.
You can choose to delete the files as follows:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"],
        path="path/to/my/target",
        clear_download=True
    )


.. |travis| image:: https://travis-ci.org/LucaCappelletti94/epigenomic_dataset.png
   :target: https://travis-ci.org/LucaCappelletti94/epigenomic_dataset
   :alt: Travis CI build

.. |sonar_quality| image:: https://sonarcloud.io/api/project_badges/measure?project=LucaCappelletti94_epigenomic_dataset&metric=alert_status
    :target: https://sonarcloud.io/dashboard/index/LucaCappelletti94_epigenomic_dataset
    :alt: SonarCloud Quality

.. |sonar_maintainability| image:: https://sonarcloud.io/api/project_badges/measure?project=LucaCappelletti94_epigenomic_dataset&metric=sqale_rating
    :target: https://sonarcloud.io/dashboard/index/LucaCappelletti94_epigenomic_dataset
    :alt: SonarCloud Maintainability

.. |sonar_coverage| image:: https://sonarcloud.io/api/project_badges/measure?project=LucaCappelletti94_epigenomic_dataset&metric=coverage
    :target: https://sonarcloud.io/dashboard/index/LucaCappelletti94_epigenomic_dataset
    :alt: SonarCloud Coverage

.. |coveralls| image:: https://coveralls.io/repos/github/LucaCappelletti94/epigenomic_dataset/badge.svg?branch=master
    :target: https://coveralls.io/github/LucaCappelletti94/epigenomic_dataset?branch=master
    :alt: Coveralls Coverage

.. |pip| image:: https://badge.fury.io/py/epigenomic-dataset.svg
    :target: https://badge.fury.io/py/epigenomic-dataset
    :alt: Pypi project

.. |downloads| image:: https://pepy.tech/badge/epigenomic-dataset
    :target: https://pepy.tech/badge/epigenomic-dataset
    :alt: Pypi total project downloads

.. |codacy| image:: https://api.codacy.com/project/badge/Grade/85bc1e3d96bf4c43a2ca70ca233a1bca
    :target: https://www.codacy.com/manual/LucaCappelletti94/epigenomic_dataset?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LucaCappelletti94/epigenomic_dataset&amp;utm_campaign=Badge_Grade
    :alt: Codacy Maintainability

.. |code_climate_maintainability| image:: https://api.codeclimate.com/v1/badges/64bfb8eb5a73959ea0d3/maintainability
    :target: https://codeclimate.com/github/LucaCappelletti94/epigenomic_dataset/maintainability
    :alt: Maintainability

.. |code_climate_coverage| image:: https://api.codeclimate.com/v1/badges/64bfb8eb5a73959ea0d3/test_coverage
    :target: https://codeclimate.com/github/LucaCappelletti94/epigenomic_dataset/test_coverage
    :alt: Code Climate Coverate
