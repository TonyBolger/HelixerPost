# HelixerPost
## Dependencies
### hd5 libraries
On fedora & co. you will need `hdf5-devel`, while on ubuntu & co. you will need `libhdf5-dev`.
Currently it's also necessary to manually build the 'lzf' compression support from h5py (https://pypi.org/project/h5py/) as a shared library and install to the hdf5 plugins directory

`tar -xzvf h5py-3.2.1.tar.gz`

`cd h5py-3.2.1/lzf/`

### Fedora
`gcc -O2 -fPIC -shared -Ilzf lzf/*.c lzf_filter.c -lhdf5 -o liblzf_filter.so`

`sudo mkdir -p /usr/local/hdf5/lib/plugin`

`sudo cp liblzf_filter.so /usr/local/hdf5/lib/plugin`

### Ubuntu
`gcc -O2 -fPIC -shared -Ilzf -I/usr/include/hdf5/serial/ lzf/*.c lzf_filter.c -lhdf5 -L/lib/x86_64-linux-gnu/hdf5/serial -o liblzf_filter.so`

`sudo mkdir /usr/lib/x86_64-linux-gnu/hdf5/plugins`

`sudo cp liblzf_filter.so /usr/lib/x86_64-linux-gnu/hdf5/plugins`



## Usage

Using 100bp sliding window for genic detection, 0.1 edge threshold, 0.8 peak threshold (before HMM)

60bp minimum CDS length per gene (after HMM)

### One Step:

`cargo run --release genome_data.h5 predictions.h5 100 0.1 0.8 60 output.gff`

### Two Step: 

`carbo build --release`

`./target/release/helixer_post_bin genome_data.h5 predictions.h5 100 0.1 0.8 60 output.gff`


