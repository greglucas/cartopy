# https://proj.org/install.html#compilation-and-installation-from-source-code

# Download and extract the specific version of PROJ
curl https://download.osgeo.org/proj/proj-7.2.1.tar.gz | tar xz

# configure and make
cd proj-7.2.1
./configure
make
make install