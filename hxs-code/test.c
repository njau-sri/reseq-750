/*This is only a test file,mainly to test whether my github can works well or not*/
#include<stdio.h>
int main(void)
{
	int i,j;
	i = 0;
	j = 0;
	for(i=1;i<10;i++)
	{
		for(j=1;j<=i;j++)
		{
			printf("%d*%d=%d\t", j, i, i*j);
		}
		printf("\n");
	}
	return 0;
}