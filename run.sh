./compile.sh
./build/lbm data/ball.txt
source env/bin/activate
python scripts/plotting2D.py output.txt