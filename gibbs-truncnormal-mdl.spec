Name:		gibbs-truncnormal-mdl
Version:	1.0
Release:	1%{?dist}
Summary:	Standalone executables for HANDLES Gibbs sampler

Group:		Applications/Engineering
License:	GPLv3
URL:		https://github.com/kralljr/handles/standalone
Source0:	gibbs-truncnormal-mdl-%{version}.tar.gz

BuildRequires:	gsl-devel >= 1.13
BuildRequires:	netcdf-devel >= 4.1.0
BuildRequires:	atlas-devel
BuildRequires:	gcc
BuildRequires:	make
BuildRequires:	pkgconfig
Requires:	    gsl >= 1.13
Requires:	    netcdf >= 4.1.0

%description
Multiply impute censored PM2.5 constituent concentrations concentrations using a likelihood-based method.


%prep
%setup -q


%build
make %{?_smp_mflags}


%install
rm -rf %{buildroot}
install -d %{buildroot}%{_bindir}
install -m 755 gibbs %{buildroot}%{_bindir}
install -d %{buildroot}%{_docdir}/handles-%{version}
install -m 644 README.md %{buildroot}%{_docdir}/handles-%{version}


%files
%{_bindir}/gibbs
%doc %{_docdir}/handles-%{version}/README.md


%changelog

