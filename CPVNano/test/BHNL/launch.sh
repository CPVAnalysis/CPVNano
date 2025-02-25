#!/bin/bash

sleep 2h
python nanoLauncher2.py --pl V42 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --sigcrab

sleep 1h
python nanoLauncher.py --pl V42_Bc --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --sigcrab

sleep 3h
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A1
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A2
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A3
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A4
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A5
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds A6
sleep 30m

python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds B2
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds B3
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds B4
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds B5
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds B6
sleep 30m

python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds C1
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds C2
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds C3
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds C4
sleep 30m

python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds D2
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds D3
sleep 30m
python nanoLauncher.py --pl V13_06Feb23 --tagnano 06Feb23 --tagflat 31Jul23 --dosignal --doflat --dosplitflat --data --ds D5
