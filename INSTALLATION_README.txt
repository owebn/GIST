conda create --name GIGN -c conda-forge -c pytorch -c nvidia python=3.9  pymol-open-source pytorch==1.13.1 torchvision torchaudio cudatoolkit=11.3 openbabel -y
conda activate GIGN
pip install scipy
pip install --no-index pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.11.0+cu113.html
pip install torch_geometric
python -m pip install PyYAML scipy "networkx[default]" biopython rdkit-pypi e3nn spyrmsd pandas biopandas
pip install scikit_learn==1.1.3
