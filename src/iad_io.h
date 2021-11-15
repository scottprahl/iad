

int Read_Header (FILE * fp, struct measure_type *m, int *params);

void Write_Header (struct measure_type m, struct invert_type r, int params);

int Read_Data_Line (FILE * fp, struct measure_type *m, int params);
