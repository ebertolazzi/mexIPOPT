### mexIPOPT and coinbrew
**by Enrico Bertolazzi**

To install and compile IPOPT use coinbrew `https://github.com/coin-or/coinbrew`.
After installation run

~~~
./coinbrew
~~~

select 1 at the first qustion

~~~
Please choose an action by typing 1-4.
 1. Fetch source code of a project and its dependencies.
 2. Build a project and its dependencies.
 3. Install a project and its dependencies.
 4. Help
=> 1
~~~

select IPOPT (11) as project to build

~~~
Please choose a main project to fetch/build by typing 1-18
or simply type the repository name of another project not
listed here.
 1. Osi
 2. Clp
 3. Cbc
 4. DyLP
 5. FlopC++
 6. Vol
 7. SYMPHONY
 8. Smi
 9. CoinMP
 10. Bcp
 11. Ipopt
 12. Alps
 13. BiCePS
 14. Blis
 15. Dip
 16. Bonmin
 17. Couenne
 18. Optimization Services
 19. MibS
 20. DisCO
 21. All
 22. Let me enter another project
=> 11
~~~

choose IPOPT version

~~~
You haven't specified a project version
It appears that the last 10 releases of Ipopt are
3.14.11
3.14.10
3.14.9
3.14.8
3.14.7
3.14.6
3.14.5
3.14.4
3.14.3
3.14.2
Do you want to work with the latest release? (y/n)
=> y
~~~

choose build (or any other directory you like) for build project

~~~
Please specify an install directory (can be relative or absolute).
=> build
~~~


**Author:**

	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
