#!/usr/bin/env bash
echo -n "Enter Eq. Distance: "
read EquiDis

RemoveLower.py Po2{s,b}AS{casscf,caspt2,mscaspt2,rassi,sorassi}.txt $EquiDis

