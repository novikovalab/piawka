#!/usr/bin/env mawk -f

# This script acts as `bedtools merge -d 1000`, joining BED records closer than 1kb to each other.
# This is done to reduce the number of tabix queries with numerous, small, adjacent BED regions.

NR==1 { 
   printf $1"\t"$2"\t"; lastchr=$1; lastend=$3; next 
}
{
   if ( $1!=lastchr || $2 - lastend > 1000 ) { 
     print lastend; printf $1"\t"$2"\t" 
   }
   lastchr=$1; lastend=$3
}
