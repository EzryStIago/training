data class1;
input ID Name $ Marks;
cards;
1     Rahul     45
1     Ajay      74
2     Ram       45
2     Girish    54
3     Simran    87
3     Priya     92
3     Riya      87
4     Tina      23
5     Dave      87
5     Ken       87
6     Albert    63
8     Alex      72
;
run;

PROC SORT DATA = class1;
BY ID;
RUN;

DATA class2;
SET class1;
BY ID;
IF FIRST.ID;
PROC PRINT;
RUN;
