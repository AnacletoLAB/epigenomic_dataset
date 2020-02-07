epigenomic_dataset
=========================================================================================
|travis| |sonar_quality| |sonar_maintainability|
|codacy| |code_climate_maintainability| |pip| |downloads|

Python package wrapping ENCODE epigenomic data
for a number of reference cell lines.

How do I install this package?
----------------------------------------------
As usual, just download it using pip:

.. code:: shell

    pip install epigenomic_dataset

Tests Coverage
----------------------------------------------
Since some software handling coverages sometime get slightly
different results, here's three of them:

|coveralls| |sonar_coverage| |code_climate_coverage|

Pipeline
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
        cell_lines=["HepG2", "H1]
    )

If you want to specify where to store the files use:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"],
        path="path/to/my/target"
    )

By default the downloaded bigWig files are not deleted.
You can choose to delete the files as follows:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"],
        path="path/to/my/target",
        clear_download=True
    )

Finally, you can use a custom **NaN threshold** to drop windows
that contain too many NaNs. For instance, if you want to drop
the rows that have more than 60% of NaN you can use:

.. code:: python

    from epigenomic_dataset import build

    build(
        bed_path="path/to/my/bed/file.bed",
        cell_lines=["HepG2", "H1"],
        nan_threshold=0.6
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
