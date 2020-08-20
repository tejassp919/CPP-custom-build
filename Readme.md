To have the the debugging functions (trace() and cotra()) you need to replace the default available stdc++.h file with the one provided.
After that go to your directory where you put your stdc++.h file and run following command
~~~
g++ -std=c++17 stdc++.h
~~~
This will create a "stdc++.h.gch" file....now rename the "stdc++.h" file to "istdc++.h" (now your code will use precompiled version of stdc++.h)

Note : Please ensure that the compiler version supports c++17, else replace all commands of "-std=c++17" with "-std=c++14"
