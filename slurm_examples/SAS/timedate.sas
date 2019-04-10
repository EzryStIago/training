/* Take date field as character string */
DATA Employee_Info;
input Emp_ID Emp_Name$ Emp_Vertical$ DOJ$;
datalines;
101 Mak SQL 18/08/2013
102 Rama SAS 25/06/2015
103 Priya Java 21/02/2010
104 Karthik Excel 19/05/2007
105 Mandeep SAS 11/09/2016
;
Run;
PROC PRINT DATA=Employee_Info;
TITLE 'USE CHARACTER STRING FOR DATE ';
Run;

/* Use INFORMAT to store date as days from Jan 1, 1960 */
DATA Employee_Info;
input Emp_ID Emp_Name$ Emp_Vertical$ DOJ;
INFORMAT DOJ ddmmyy10.;
datalines;
101 Mak SQL 18/08/2013
102 Rama SAS 25/06/2015
103 Priya Java 21/02/2010
104 Karthik Excel 19/05/2007
105 Mandeep SAS 11/09/2016
;
Run;
PROC PRINT DATA=Employee_Info;
TITLE 'USE INFORMAT';
Run;

/* Use INFORMAT to store date as days from Jan 1, 1960  */
/* Use FORMAT to display data fields in ddmmyyyy format */
DATA Employee_Info;
input Emp_ID Emp_Name$ Emp_Vertical$ DOJ;
INFORMAT DOJ ddmmyy10.;
FORMAT DOJ ddmmyy10.;
datalines;
101 Mak SQL 18/08/2013
102 Rama SAS 25/06/2015
103 Priya Java 21/02/2010
104 Karthik Excel 19/05/2007
105 Mandeep SAS 11/09/2016
;
Run;
PROC PRINT DATA=Employee_Info;
TITLE 'USE BOTH FORMAT AND INFORMAT';
Run;
