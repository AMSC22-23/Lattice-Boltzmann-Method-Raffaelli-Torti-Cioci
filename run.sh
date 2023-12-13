./compile.sh
./build/lbm data/lid-driven-cavity.txt 0.002 200
source env/bin/activate
python scripts/plotting2D.py output.txt