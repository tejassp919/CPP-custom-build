template <typename T>               void __p (T a);
template <typename T, typename F>   void __p (pair<T, F> a);
template <typename S, size_t T>     void __p (array<S, T> arr);
template <typename T>               void __p (vector<T> &a);

template <typename T>
void __p (T a) { cout << a << " "; }

template <typename T, typename F>
void __p (pair<T, F> a) { 
    cout << "{ "; 
    __p(a.first); __p(a.second); 
    cout << "} "; 
}

template <typename S, size_t T>
void __p (array<S,T> arr) {
    for(int i = 0; i < T; i++) cout << arr[i] << " ";
    cout << endl;
}

template <typename T>
void __p (vector<T> &a) { 
    cout << endl << " { "; 
    for (auto p : a) __p(p); 
    cout << "}"; 
}

template <typename Arg1>
void __t(const char *name, Arg1 &&arg1) { 
    cout << name << " : { "; 
    for (auto p : arg1) __p(p); 
    cout << "}" << endl; 
}

template<typename Arg1, typename ... Args>
void __t(const char *names, Arg1 &&arg1, Args &&... args) {
    int bracket = 0, i = 0;
    for(; ; i++)
        if (names[i] == ',' && bracket == 0)
            break;
        else if (names[i] == '(')
            bracket++;
        else if (names[i] == ')')
            bracket--;
    cout.write(names, i) << " : ";
    __p(arg1);  cout << endl;
    __t(names + i + 1, args...);
}

template <typename Arg1>
void __f(const char *name, Arg1 &&arg1) {
    cout << name << " : "; 
    __p(arg1); 
    cout<<endl;
}

template <typename Arg1, typename ... Args>
void __f(const char *names, Arg1 &&arg1, Args &&... args) {
    int bracket = 0, i = 0;
    for(; ; i++)
        if (names[i] == ',' && bracket == 0)
            break;
        else if (names[i] == '(')
            bracket++;
        else if (names[i] == ')')
            bracket--;
    cout.write(names, i) << " : ";
    __p(arg1);  cout << "| ";
    __f(names + i + 1, args...);
}

template <typename Arg1, typename Arg2>
void __f(const char *names, Arg1 arg1[], Arg2 &&arg2){
    int i = 0;
    for(; ; i++) if (names[i] == ',') break;
    cout.write(names, i) << " : { ";
    for(int i = 0; i < arg2; i++) __p(arg1[i]);
    cout << "} " << endl;
}

#define cotra(...) { cout<<"Line:"<<__LINE__<<" | "; __t(#__VA_ARGS__, __VA_ARGS__); }
#define trace(...) { cout<<"Line:"<<__LINE__<<" | "; __f(#__VA_ARGS__, __VA_ARGS__); }
int begtime = clock();
#define end_routine() cout << "\n\nTime elapsed: " << (clock() - begtime)*1000/CLOCKS_PER_SEC << " ms\n\n";
