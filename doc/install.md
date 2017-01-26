# Install

CosmoSlik has the following requirements:

* Python 3
* Cython >= 0.16
* numpy >= 1.11.0
* scipy >= 0.17.0
* setuptools >= 3.2

If you don't have these, they will be installed along with CosmoSlik, although its probably best to try and install them beforehand with whatever package manager you use on your system. 

To install CosmoSlik, first get the source code<sup>[1](#fn)</sup>, 

```bash
git clone https://github.com/marius311/cosmoslik.git
cd cosmoslik
```

then run, 

```bash
python setup.py develop --user
```

This installs CosmoSlik in-place so you can update it later by simlpy running `git pull`.

<a name="fn">[1]</a>: If your system doesn't have `git` installed, you can also manually download the source code from [https://github.com/marius311/cosmoslik](https://github.com/marius311/cosmoslik)
