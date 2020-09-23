Parallel Processing in Sidle
----------------------------

Sidle is built to be run in parallel using `dask`_. This sysstem allows for flexible implementations. All side commands except ``reconstruct taxonomy`` allow parallel processing. There are three ways to manage parallel processing:

* ``--debug`` which turns off all parallel processing. This is primarily implemented for testing
* ``--p-n-workers`` will create a new dask cluster with the specified number of workers and will use the avalaible resources
* ``--p-client-address`` allows you to pass in a pre-configured cluster for processing. You can learn more about setting up `dask clients in their documentation`_.

This hopefully provides a flexible interface for parallel procesisng of your samples

.. websites
.. _dask: https://dask.org/
.. _dask clients in their documentation: https://docs.dask.org/en/latest/setup.html#