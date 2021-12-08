; Atom type definitions for use with CanvasMCS

; Atoms with the same atomic number are equivalent.
; All bonds are equivalent.
; Distinguish some carbons.

; Application of atom typing rules, least specific to most specific.

; Current RMSD calculation only works if all classes have same atomic number.

[*]         > 1  ;
[#1]        > 2  ;
[#6]        > 3  ;
[#6;r5;Cx4] > 4  ; 
[#6;r6]     > 4  ;
c1ccccc1    > 5  ;
[CR0]       > 6  ;
[#7]        > 7  ;
[#7;r5]     > 12 ;
[#8]        > 13 ;
O=*         > 14 ;
[#8;r5]     > 15 ;
[#15]       > 17 ;
[#16]       > 18 ;
[#16;r5]    > 19 ;
[#9]        > 20 ;
[#17]       > 20 ;
[#35]       > 20 ;
[#53]       > 20 ;

