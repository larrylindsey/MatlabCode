//For testing image processing utilities library

#include <stdio.h>

#include <file_input_output.h>

int main(int argc, char * argv[]){
  printf("\nSTART: testbed for image processing utilities\n");
  {
    char test_file[] = "test.raw";
    iput::io::test_file_input_output(test_file);
  }
  printf("\nSTOP: testbed for image processing utilities\n");
}
