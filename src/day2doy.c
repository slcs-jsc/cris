#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  int doy;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <year> <mon> <day>");

  /* Read arguments... */
  const int year = atoi(argv[1]);
  const int mon = atoi(argv[2]);
  const int day = atoi(argv[3]);

  /* Convert... */
  day2doy(year, mon, day, &doy);
  printf("%d %d\n", year, doy);

  return EXIT_SUCCESS;
}
