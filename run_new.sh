g++ src/Fibo.cpp -o executable -Ofast
./executable
source env/bin/activate
python scripts/plotting2D.py output.txt