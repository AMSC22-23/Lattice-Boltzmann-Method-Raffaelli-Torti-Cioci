./compile.sh
./build/lbm data/lid-driven-cavity.txt
source env/bin/activate
python scripts/plotting2D.py output.txt