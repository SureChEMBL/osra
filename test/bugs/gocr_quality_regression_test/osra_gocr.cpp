
#include <iostream>

#include <stdlib.h>
#include <string.h>

extern "C" {
#include <pgm2asc.h>
}

using namespace std;

/* Actual max height is 12, but we leave some more for extensions: */
const char* TESTS[][50] =
{
  /* These show where 0.45 is better */
  /* Test1: "3" is not detected */
  {
    "##############",
    "####......####",
    "###........###",
    "###..#......##",
    "########...###",
    "########...###",
    "########...###",
    "########..####",
    "######....####",
    "#####.....####",
    "#######.....##",
    "#########....#",
    "#########....#",
    "##########...#",
    "#########...##",
    "##...###....##",
    "#..........###",
    "#.........####",
    "###..#.#######"
  },
  /* Test2: "3" is not detected */
  {
    "##############",
    "####......####",
    "###........###",
    "###..#......##",
    "########...###",
    "########...###",
    "########...###",
    "########..####",
    "######....####",
    "#####.....####",
    "#######.....##",
    "#########....#",
    "#########....#",
    "##########...#",
    "#########...##",
    "##...###....##",
    "#..........###",
    "#.........####",
    "###..#.#######",
    "##############",
  },
  /* Test3: "3" is not detected */
  {
    "############",
    "###......###",
    "###.......##",
    "###........#",
    "#######....#",
    "#######....#",
    "#######...##",
    "#######...##",
    "######....##",
    "#####.....##",
    "######.....#",
    "#######.....",
    "########....",
    "########...#",
    "########...#",
    "#...###....#",
    "#.........##",
    "#........###",
    "###.....####",
    "############"
  },
  /* The rest show where 0.48 is better */
  /* Test4: nothing should be detected */
  {
    "#######################",
    "#######################",
    "#######################",
    "#########..#..#########",
    "########.........######",
    "########........#######",
    "#########......########",
    "##########....#########",
    "###########...#########",
    ".##########...#########",
    ".##########...#########",
    ".#########....#########",
    "...#..........#########",
    "..............#########",
    "..............#########",
    ".########.....#########",
    ".#########....#########",
    ".#########....#########",
    "##########....#########",
    "##########....#########",
    "#########.....#########",
    "########......#########",
    "#######........########"
  },
  /* Test5: "N" should be detected */
  {
    "###############",
    "#....#####....#",
    "##....#####....",
    "##.....####..##",
    "##.....#####.##",
    "##...#..###..##",
    "##..##...##..##",
    "##..###...#..##",
    "##..####.....##",
    "##..#####....##",
    "##..#####....##",
    "##..######...##",
    "#.....#####..##",
    "#....######..##",
    "###############"
  },
  /* Test6: "C" should be detected */
  {
    "########.............",
    "#####..#.............",
    "#####................",
    "#####................",
    "####.................",
    "####........#########",
    "####........#########",
    "##..........#########",
    "##..........#########",
    "##.........##########",
    "##.........##########",
    "##.........##########",
    "##.........##########",
    "##.........##########",
    "##.........##########",
    "##.........##########",
    "##..........#########",
    "##.#........#########",
    "#.#.........#########",
    "###............##.##.",
    "####............#..#.",
    "####.................",
    "#####................"
  }
};

job_t *JOB;
job_t *OCR_JOB;

char run_test(int n)
{
  int height = 0;
  int width = strlen(TESTS[n][0]);

  while (TESTS[n][height] != NULL)
    {
      height++;
    }

  const char** image = TESTS[n];

  cout << "Test " << n + 1 << ": width x height = " << width << "x" << height << endl;

  job_t job;

  job_init(&job);
  job_init_image(&job);

  job.cfg.cfilter = (char *) "oOcCnNHFsSBuUgMeEXYZRPp23456789";
  job.src.p.x = width;
  job.src.p.y = height;
  job.src.p.bpp = 1;
  job.src.p.p = (unsigned char *) malloc(job.src.p.x * job.src.p.y);

  for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col++)
        job.src.p.p[row * width + col] = image[row][col] == '#' ? 255 : 0;
    }

  JOB = &job;
  OCR_JOB = &job;

  try
    {
      pgm2asc(&job);
    }
  catch (...)
    {
    }

  char *l = (char *) job.res.linelist.start.next->data;

  char c = 0;

  if (l != NULL)
    c = l[0];

  if (isalnum(c))
    {
      // Character recognition succeeded for GOCR:
      cout << "Found c=" << c << endl;
    }
  else
    {
      cout << "Failed c=" << c << endl;
    }
}

int main()
{
  for (unsigned int n = 0; n < 6; n++)
    {
      run_test(n);
    }
}
