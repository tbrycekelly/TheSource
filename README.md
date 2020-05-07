# TheSource

__TheSource__ is an R package for the rest of us, whose sole mission is to allow us to spend more time doing what we want to do and less time figuring out how to get there. Included are numerous helper functions, such as a function to convert between excel datetime objects and "real" datetime objects. In addition, a number of oceanographic functions are included, which I have used extensively in my research and in others. Examples of oceanographic functions include the calculation of seawater density from T and S, or calculating air-sea gas exchange coefficeints.

Don't worry of any of these functions don't sound useful to you, there is plenty of good stuff to go around. So grab a beer and check out some of the tutorials!

_To install:_
```R
install.packages('devtools')
devtools::install_github('tbrycekelly/TheSource')
```

To update an existing installation of TheSource, simply run:
```R
TheSource.update()
```

A few recommendations for Windows and Mac users:

1. On Windows/PC, it's a good idea to install __Rtools__ from the R website. This will allow you to compile some packages from source.
2. For Mac, the _commandline tools_ will be helpful for the same reason. In Terminal run __xcode-select --install__.


## Citation

If you use the library in a publication please include the following citation:


Kelly, T.B., _TheSource: It will hold your beer. Zenodo (2019). DOI: 10.5281/zenodo.3468524
[![DOI](https://zenodo.org/badge/209631718.svg)](https://zenodo.org/badge/latestdoi/209631718)


Cheers,
T

<img src="https://github.com/tbrycekelly/TheSource/blob/master/logo.png" width="250">
