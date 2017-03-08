#include "Main.h"

void main(void)
{
	repeater();
	system("pause");
}
void starter(int Eqn)
{
	switch (Eqn)
	{
	case 1:
		printf("Parabolic\n");
		Parabolic::main();
		break;
	case 2:
		printf("Elliptic\n");
		Elliptic::main();
		break;
	case 3:
		printf("Hyperbolic\n");
		Hyperbolic::main();
		break;
	}
}
int Eqnselector(void)
{
	int Eqn;
	printf("Parabolic: 1\n");
	printf("Elliptic: 2\n");
	printf("Hyperbolic: 3\n");
	while (1)
	{
		printf("What equation you want to solve?: ");
		scanf("%d", &Eqn);
		getchar();
		if (Eqn >= 1 && Eqn <= 3)
		{
			break;
		}
		else
			printf("Input proper value\n");
	}
	printf("\n");
	return Eqn;
}
void repeater(void)
{
	int checker=1;
	do{
		int Eqn = Eqnselector();
		starter(Eqn);
		printf("\nOnce more?\n");
		printf("[yes : 1, no : 0]: ");
		scanf("%d", &checker);
		getchar();
		if (checker != 1 && checker != 0)
		{
			printf("Wrong Input, program break\n");
			return;
		}
		printf("\n");
	}while (checker);
}