# capi

`capi` is a simple Go package that can be compiled as a C shared library.
The `Makefile` shows how it should be compiled.

`test-root.C` shows how it can be used from C/C++ ROOT.

## example

```
$> make clean
/bin/rm -f ./libcapi.so ./libcapi.h

$> make
go build -v -o libcapi.so -buildmode=c-shared
runtime/cgo
github.com/rmadar/go-topquark-reco/capi

$> make test
Processing ./test-root.C...
entries: 100
entry: 0
entry: 10
entry: 20
evt: 195742678 (entry=25)
top[0]: (-108.192, -6.24721, -168.478, 264.359)
top[1]: (6.55442, 65.1181, 171.649, 251.998)
entry: 30
entry: 40
evt: 195742507 (entry=40)
top[0]: (88.8442, 39.3959, 146.573, 246.344)
top[1]: (-96.9154, -60.9951, 16.3345, 207.692)
entry: 50
evt: 195741255 (entry=51)
top[0]: (83.4215, -75.7093, 70.2792, 217.684)
top[1]: (-70.996, 52.4966, 13.1201, 194.228)
evt: 195741196 (entry=55)
top[0]: (-36.2465, -158.237, -226.089, 327.453)
top[1]: (40.3798, 155.862, -163.054, 286.821)
entry: 60
entry: 70
evt: 195742754 (entry=72)
top[0]: (28.5613, 125.463, 200.766, 294.312)
top[1]: (-39.6241, -111.878, 287.533, 355.694)
entry: 80
entry: 90
evt: 195742040 (entry=92)
top[0]: (25.9241, 38.6864, -21.4179, 179.955)
top[1]: (4.34166, 24.7786, 231.494, 289.791)
evt: 195741784 (entry=95)
top[0]: (-35.6923, -32.4248, -42.3103, 184.043)
top[1]: (55.0593, 17.5964, 173.504, 251.398)
n-reco: 7
```
