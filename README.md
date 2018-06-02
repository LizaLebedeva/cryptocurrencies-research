# cryptocurrencies-research
This is part of my master thesis on topic *"Spillovers between cryptocurrencies. Network map of cryptocurrencies"*.

The aim of thesis is to determine the interaction between different coins and construct network map of cryptocurrencies. To do this, I applied methodology of Variance decomposition based on VAR time series analysis.

To analyze connectedness of cryptocurrencies market:
* Define coins for analysis (your own list or using [get-top-coins.ipynb](https://github.com/LizaLebedeva/cryptocurrencies-research/blob/master/get-top-coins.ipynb) )
* Get historical data ([crawl-data.ipynb](https://github.com/LizaLebedeva/cryptocurrencies-research/blob/master/crawl-data.ipynb))
* Run [return-network-analysis.R](https://github.com/LizaLebedeva/cryptocurrencies-research/blob/master/return-network-analysis.R) (be sure that you also have [functions.R](https://github.com/LizaLebedeva/cryptocurrencies-research/blob/master/functions.R)) which provides with different total connectedness measurements of a network (static, dynamic), pairwise connectedness of coins, most connected and less connected coins, individual analysis of particular coin.
