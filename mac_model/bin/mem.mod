V34 :0x24 mem
42 /home/gleung/ats730/mac_model//src/mem.f90 S624 0
01/29/2024  17:17:09
use model_vars public 0 direct
use grid_constants private
enduse
S 624 24 0 0 0 6 1 0 5013 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 mem
S 626 23 0 0 0 6 640 624 5032 4 0 A 0 0 0 0 B 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 nz
R 640 6 1 grid_constants nz
S 793 23 5 0 0 0 794 624 6176 0 0 A 0 0 0 0 B 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 allocate_mem
S 794 14 5 0 0 0 1 793 6176 0 400000 A 0 0 0 0 B 0 9 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 9 0 624 0 0 0 0 allocate_mem
F 794 0
S 795 23 5 0 0 0 796 624 6189 0 0 A 0 0 0 0 B 0 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 deallocate_mem
S 796 14 5 0 0 0 1 795 6189 0 400000 A 0 0 0 0 B 0 33 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 33 0 624 0 0 0 0 deallocate_mem
F 796 0
Z
Z
