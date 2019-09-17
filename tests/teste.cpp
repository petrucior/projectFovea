#include<iostream> 
#include<vector>
using namespace std; 

// An abstract class with constructor 
class Base { 
protected: 
  vector<int> x, y;
public: 
  virtual void fun() = 0;
  Base(vector<int> i, vector<int> j) { x = i; y = j; } 
};

class Derived: public Base{
public: 
  Derived(vector<int> a, vector<int> b);
  void fun() { cout << "x = " << x[0] << ", p = " << y[0] << endl; } 
};

Derived::Derived( vector<int> a, vector<int> b ) : Base( a, b ){}

int main(void){
  vector<int> v1, v2;
  v1.push_back( 4 );
  v2.push_back( 5 );
  Derived d(v1, v2); 
  d.fun(); 
  return 0; 
} 
