#
# spec file for package OSRA
#
# Copyright (c) 2009 SUSE LINUX Products GmbH, Nuernberg, Germany.
#
# All modifications and additions to the file contributed by third parties
# remain the property of their copyright owners, unless otherwise agreed
# upon. The license for this file, and modifications and additions to the
# file, is the same license as for the pristine package itself (unless the
# license for the pristine package is not an Open Source License, in which
# case the license is the MIT License). An "Open Source License" is a
# license that conforms to the Open Source Definition (Version 1.9)
# published by the Open Source Initiative.

# Please submit bugfixes or comments via http://bugs.opensuse.org/
#

# norootforbuild

# %{name} and %{version} should be either defined here or inherited from ~/.rpmmacros

Name:			%{name}
BuildRequires:	glibc-devel, libstdc++43-devel, tclap >= 1.2, potrace-devel >= 1.8, gocr-devel >= 0.49, ocrad-devel >= 0.20, libopenbabel-devel >= 2.2, libGraphicsMagick++-devel >= 1.2.5, docbook-xsl-stylesheets => 1.74.0, libxslt
Url:			http://osra.sourceforge.net/
Summary:		A command line chemical structure recognition tool
Version:		%{version}
Release:		1.0
Group:			Productivity/Graphics/Other
Requires:		potrace-lib >= 1.8, libopenbabel3 >= 2.2, libGraphicsMagick++2 >= 1.2.5
#Provides:
License:		GPL v2 or later
Source0:		%{name}-%{version}.tar.gz
#Patch0:		Makefile.in.patch
BuildRoot:		%{_tmppath}/%{name}-%{version}-build

%description
OSRA is a utility designed to convert graphical representations of chemical structures into SMILES or SDF.
OSRA can read a document in any of the over 90 graphical formats parseable by GraphicMagick and generate
the SMILES or SDF representation of the molecular structure images encountered within that document.

Authors:
--------
    Igor Filippov <igorf@helix.nih.gov>


%prep
%setup -n %{name}-%{version}
#%patch0 -p0

%build
./configure --prefix=/usr --enable-docs --datadir='${datarootdir}/${PACKAGE_NAME}' --docdir='${datarootdir}/doc/packages/${PACKAGE_NAME}'
make CXXFLAGS='${RPM_OPT_FLAGS}'

%install
make install DESTDIR=$RPM_BUILD_ROOT

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-, root, root)
/usr/bin/osra
/usr/share/osra
/usr/share/man/man1/osra.1*
/usr/share/doc/packages/osra

%doc README

# spec file ends here

%changelog
* Thu Jul 11 2010 dma_k@mail.ru
- Initial SuSE package
