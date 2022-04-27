
from dask.distributed import Client


def _setup_dask_client(debug=False, cluster_config=None, n_workers=1,
    address=None):
    """
    Sets up a Dask client and daskboard

    Parameters
    ----------
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    client_config: dict, optional
        A dictionary describing configuration parameters for the dask client.
        More information about configuring the dask scheduler and dask client 
        can be found at
            https://docs.dask.org/en/latest/setup/single-distributed.html
        The client_config sueprceeds the n_workers value, so if you want 
        multi threading, that should be specified here.
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avalaibel resources.
    address: str, optional
        The IP address for the client
    """

    if debug:
        pass

    elif cluster_config is not None:
        client = Client(**client_config.to_dict())

    elif address is not None:
        client = Client(address)
        
    else:
        client = Client(n_workers=n_workers, processes=True)
