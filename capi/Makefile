.PHONY: all build

GO := go1.15.6
RM := /bin/rm

build: libtopreco.so

clean:
	$(RM) -f ./libtopreco.so ./libtopreco.h

libtopreco.so: topreco.go
	$(GO) build -v -o libtopreco.so -buildmode=c-shared

test: build test-root.C
	root.exe -b -q ./test-root.C
