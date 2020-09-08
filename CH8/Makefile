R = R CMD SHLIB
CC = gcc
SFX1 = .so
SFX2 = .o

.PHONY: all

ALL : $(addsuffix $(SFX1),KnS1998)

$(addsuffix $(SFX1),KnS1998) : KnS1998.c $(addsuffix $(SFX2),cubic_spline) $(addsuffix $(SFX2),linear_interpolation) $(addsuffix $(SFX2),which)
	$(R) $^

$(addsuffix $(SFX2),cubic_spline): cubic_spline.c
	$(CC) -c $^

$(addsuffix $(SFX2),linear_interpolation): linear_interpolation.c
	$(CC) -c $^

$(addsuffix $(SFX2),which): which.c
	$(CC) -c $^

.PHONY: clean

clean:
	rm $(addsuffix $(SFX1),KnS1998) $(addsuffix $(SFX2),KnS1998) $(addsuffix $(SFX2),cubic_spline) $(addsuffix $(SFX2),which) $(addsuffix $(SFX2),linear_interpolation) 

