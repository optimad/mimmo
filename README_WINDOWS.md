#  mimmo on Microsoft Windows

This is a short guide to set a mimmo compliant 64bit Windows environment using MSYS2/MingGW64.
This guide is linked to bitpit related one (bitpit README_WINDOWS.md): all basic steps
(also needed by mimmo) are summed up there.
If a successfull bitpit installation on Windows OS was already performed, mimmo installation is almost
straighforward.

## Requirements:
- <B>bitpit</B> (compliant Windows bitpit installation and its dependencies)

- <B>MinGW64 packages</B>:

  - mingw-w64-x86_64-cgns (optional, needed for iocgns module activation)
  - mingw-w64-x86_64-parmetis (parMETIS and METIS installation, needed by MPI version only)

## Procedure
#### Install mimmo dependencies

**__bitpit__**
Follow the guide README_WINDOWS.md coming with bitpit distribution.  

**__CGNS__**
If the *iocgns* module of mimmo is required, open a *mingw64* shell.

Install *cgns* package  with:
```bash
user@machine MINGW64 ~
> pacman -S mingw-w64-x86_64-cgns
```

**__parMETIS and METIS__**
The libraries are needed only in MPI version of mimmo. For serial compiling skip this part.

If not, open a *mingw64* shell. Install *parmetis* and its sub-dependency *metis* with:
```bash
user@machine MINGW64 ~
> pacman -S mingw-w64-x86_64-parmetis
```

#### mimmo
- Download *mimmo* master archive at https://github.com/optimad/mimmo/archive/master.zip or git clone it using SSH or HTTPS.

- Use the *MinGW64* development environment to install __*mimmo*__, following the  installation instructions
available at https://github.com/optimad/mimmo/blob/master/INSTALL.md or locally in the *INSTALL.md* file in *mimmo* root folder.

- In order to compile a Windows native version of *mimmo* you have to:

  - specify the __*MinGW Makefiles Generator*__ to cmake, e.g. in your mimmo build folder
  ```bash
  user@machine MINGW64 ~
  > cmake -G "MinGW Makefiles" ..
  ```
  or if *cmake-gui* is used, enable a __*MinGW Makefiles*__ project.

  - use __*mingw32-make*__ in place of standard *make* command to build *mimmo*.


- **__Warning__**: Despite mimmo's Linux version supports IOOFOAM module, for I/O handling
of OpenFoam meshes, no module is available in the Windows twin version. In our knowledge,
no official OpenFoam distro running on Windows expose development libraries in the format
needed for a direct MIngW compiling/linking in a third code. However, any idea, suggestion,
workaround or share-of-experience is most welcome on the matter. We will be glad to use your tips
to enrich the current guide.   
