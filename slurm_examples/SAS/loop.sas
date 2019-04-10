data A; 
do i = 1 to 4;   
y = i**2; /* values are 2, 5, 9, 16, 25 */  
output; 
end; 
run;

PROC PRINT DATA=A;
	TITLE 'do loop number 1';
RUN;


data B;
do i = 1 to 5 by 0.3;  
y = i**2; /* values are 1, 2.25, 4, ..., 16, 20.25, 25 */ 
output;
end; 
run;

PROC PRINT DATA=B;
	TITLE 'do loop number 2';
RUN;

