#!/usr/bin/env ruby
#---------------update
require '~/.scripts/pool.rb'
def f(i)
	system("echo #{i}")
	system("g++ -std=c++11 nlopt_turm_mwmix_232_probability.cc matrix.cc -O2 -o 3.out ")
	system("cat work.sh |\
			 sed \"/^i=/c \i=#{i}\ \" >  work_#{i}.sh")
	system("chmod +x work_#{i}.sh")
	system("./work_#{i}.sh ")
	system("rm -rf work_#{i}.sh ")
end
#---------------------------------------------------------------
pool = Pool.new(20,5)
N =20
(1..N).each do |i|
   	pool.schedule {f(i)}
end
pool.shutdown
sleep 2
p "All finished, congratulations"
