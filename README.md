# Standalone Gibbs Sampler

The standalone Gibbs sample implements a likelihood-based method for imputing values below the detection limit. The method assumes the complete, unobserved data are multivariate normal and uses an MCMC approach to estimate the parameters. We use conjugate priors for the mean and covariance of the multivariate normal distribution and directly sample from the full conditional distributions using a Gibbs sampler. The full conditional for the censored data is truncated normal, where the data are truncated above by the MDL.

See <https://github.com/kralljr/handles> for additional details.

## Installation

The standalone Gibbs sampler requires two external libraries: [GNU Scientific Library (NetCDF)](http://www.gnu.org/software/gsl/) for numerical functions and [Network Common Data Form (NetCDF)](http://www.unidata.ucar.edu/software/netcdf/) for storage. Precompiled binaries exist for many operating systems, but if absent, the libraries may be compiled from source.

When building the Gibbs sampler from source, [pkg-config](http://www.freedesktop.org/wiki/Software/pkg-config/) is used by default to manage compiler flags, but is not required. A working C compiler and `make` are also required to build from source.

Once the dependencies are installed, to build the Gibbs sampler, review the options in `Makefile`, then invoke `make`. The compiler flags have been tested with `pkg-config` present, GSL and NetCDF installed into a standard system library (`/usr` or `/usr/local`), and with GCC as the compiler.

## Usage

Once built, the Gibbs sampler may be used:

    ./gibbs DATA MDLS OUTPUT

where DATA is the name of an input data file, MDLS is a the name of the file containing MDLS, and the OUTPUT is the name of the output NetCDF file. To run the executable from somewhere other than the current directory, replace `./` with the name of the directory or install the executable to a directory on the users' path and omit the leading `./`.

The executable takes the optional flags:

    -n int   number of iterations (default 1000)
    -b int   number of iterations to burn-in (default 0)
    -d int   number of random draws (default 1)
    -r int   random seed (default 0)
    -h       show this help

## Data

The input data (PM concentations) and minimum detection limits (MDLs) are CSV files. Each column of the data (or MDLs) represents a particulate and each row contains the concentration (or MDL) for that particulate on a given day. Column labels may be provided, but row labels must be absent.

## Operating system specific installation notes

### OS X

GSL and NetCDF are available in both [Homebrew](http://brew.sh/) and [MacPorts](https://www.macports.org/).

Edit `Makefile` to enable the Accelerate framework within the linker flags (`LDFLAGS`).

### Red Hat Enterprise Linux 6

GSL is available as an RPM from the RHEL base repository. NetCDF may be built from source or is available from third party repositories such as [RepoForge](http://pkgs.repoforge.org/netcdf/).

### Fedora

Both GSL and NetCDF are available from the Fedora base repository.

### Installing GSL from source

Directions are an overview only. Please refer to the file INSTALL within a source GSL package for details. Directions are given for GSL 1.16, but previous versions may work as well.

Download GSL from <http://ftpmirror.gnu.org/gsl/gsl-1.16.tar.gz>.

Move the tarball to the current directory or change the current directory to the location of the downloaded tarball. Unpack the tarball:

    tar xzf gsl-1.16.tar.gz

Change directory to the unpacked tarball:

    cd gsl-1.16

Configure, build, and install the package:

    ./configure
    make
    make install

By default, the package is installed in `/usr/local`. To change this location to `[SOME DIRECTORY]`, pass the option `--prefix=[SOME DIRECTORY]` to the configuration script:

    ./configure --prefix=[SOME DIRECTORY]

then build and install as before.

    make
    make install

### Installing NetCDF from source

Directions are an overview only. The resulting NetCDF library is sufficient to use the Gibbs executable, but may not be usable for other applications requiring NetCDF. Please see <http://www.unidata.ucar.edu/software/netcdf/docs/building.html> for additional details.

Download NetCDF from <ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.2.tar.gz>.

Move the tarball to the current directory or change the current directory to the location of the downloaded tarball. Unpack the tarball:

    tar xzf netcdf-4.3.2.tar.gz

Change directory to the unpacked tarball:

    cd netcdf-4.3.2

Configure, build, and install the package:

    ./configure --disable-dap
    make check install

By default, the package is installed in `/usr/local`. To change this location to `[SOME DIRECTORY]`, pass the option `--prefix=[SOME DIRECTORY]` to the configuration script:

    ./configure --prefix=[SOME DIRECTORY] --disable-dap

then build and install as before.

    make check install

## Copying

The standalone Gibbs sampler is copyright 2012-2014 and licensed under the GPLv3.
