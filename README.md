# faultToVect
> A C++ application to transform faults raster prediction to vector format.

[![CircleCI][circleci-badge]][circleci-url]
[![CodeFactor Grade][codefactor-badge]][codefactor-url]
[![Documentation][documentation-badge]][documentation-url]
[![License][license-badge]][license-url]

## Getting Started

To build the project for windows:
- Install qgis from [osgeo4w][osgeo4w-url]
- Install vcpkg to help install remaining dependencies [vcpkg][vcpkg-url]
- Set following environment variables in powershell script file (set_env.ps1): 
	- OSGEO4W_ROOT
	- VS17COMNTOOLS ( probably "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\Common7\Tools" )
	- VCPKG_ROOT
- Install dependencies by running `path\to\vcpkg.exe --triplet x64-windows install "@.vcpkg_deps.txt"`
- Run cmake command in build folder `cmake .. "-DCMAKE_TOOLCHAIN_FILE=path\to\vcpkg\scripts\buildsystems\vcpkg.cmake" -G "Visual Studio 16 2019"`
- Build `~/PROJECTNAME/build/cmake --build . --config Release`
- Execute the tests `./out/build/ctest`
- You can execute the program by `./out/build/Release/faultToVect.exe`

To update the docker image:
- Edit the Dockerfile to your needs
- Build docker image `sudo docker build -t IMAGENAME .`
- Tag docker image with dockerhub username `sudo docker tag IMAGENAME:TAG DOCKERHUBUSERNAME/IMAGENAME:TAG`
- Push docker image to dockerhub `sudo docker push DOCKERHUBUSERNAME/IMAGENAME:TAG`

To change/add dependencies:
- Edit the vcpkg part of `.cirlceci/config.yml` to your needs
```
- run:
    name: Install vcpkg dependencies
    command: ./../../vcpkg/vcpkg install DEPENDENCIES
```

### Prerequisites/Dependencies

- [cmake][cmake-url] – Open-Source, cross-platform build tool
- [vcpkg][vcpkg-url] – C++ Library Manager for Windows, Linux, and MacOS
- [python 3][python-url] – A programming language used to convert ctest results with a xml transformation (xslt)


## Acknowledgments

- Converting CTest results int JUnit XML – https://stackoverflow.com/a/21688776/1541782
- Doxygen GitHub-Action – https://github.com/mattnotmitt/doxygen-action
- gh-pages GitHub-Action – https://github.com/peaceiris/actions-gh-pages
- Dockerfile Tips – https://blog.container-solutions.com/6-dockerfile-tips-official-images

[circleci-url]: https://circleci.com/gh/Ben1980/cpptemplate
[codefactor-url]: https://www.codefactor.io/repository/github/ben1980/cpptemplate
[documentation-url]: https://ben1980.github.io/cpptemplate/
[license-url]: https://github.com/Ben1980/cpptemplate/blob/master/LICENSE
[circleci-badge]: https://img.shields.io/circleci/build/gh/Ben1980/cpptemplate
[codefactor-badge]: https://img.shields.io/codefactor/grade/github/ben1980/cpptemplate
[documentation-badge]: https://img.shields.io/github/workflow/status/Ben1980/cpptemplate/Documentation?label=Documentation
[license-badge]: https://img.shields.io/github/license/Ben1980/cpptemplate
[cmake-url]: https://cmake.org/
[fmt-url]: https://fmt.dev/latest/index.html
[doctest-url]: https://github.com/onqtam/doctest
[rep-url]: https://github.com/Ben1980
[linkedin-url]: https://www.linkedin.com/in/benjamin-mahr-728a1639/
[twitter-url]: https://twitter.com/BenMahr
[mail]: ben.amhr@gmail.com
[vcpkg-url]: https://github.com/microsoft/vcpkg
[osgeo4w-url]: https://www.osgeo.org/projects/osgeo4w/
[python-url]: https://www.python.org/

[v1.0.0]: https://github.com/Ben1980/cpptemplate/releases/tag/v1.0.0
