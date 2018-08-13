# MOB-suite: Software tools for clustering, reconstruction and typing of plasmids from draft assemblies

```/Users/kirill/WORK/MOBSuiteWrapper/galaxy_tools/tools/mob_suite```

These are Galaxy wrappers for [MOB-Suite](https://github.com/phac-nml/mob-suite): a set of tools designed to be a modular set of tools for the typing and reconstruction of plasmid sequences from WGS assemblies.


For more information please refer to https://github.com/phac-nml/mob-suite
 
#Running planemo on nyx
Running on nyx server planemo from the ```/home/kbessonov/galaxy_tools/mob_suite_wrapper``` directory.
 
#Deploying docker image of planemo

```docker run --privileged=true -p 8010:80 -p 9009:9009 -v ~/WORK/MOBSuiteWrapper/galaxy_tools/tools:/opt/galaxy/tools -it --rm planemo/interactive

docker run --privileged=true -p 8010:80 -p 9090:9090 -v ~/WORK/MOBSuiteWrapper/galaxy_tools/tools:/opt/galaxy/tools -i -t bgruening/planemo

docker pull bgruening/planemo

```

Check reference [here](https://github.com/galaxyproject/planemo/blob/master/docs/appliance.rst)


Attach to instance
```
docker exec -it 2943fb75f22d /bin/bash
>root@e43a62031777:/home/ubuntu# su - ubuntu

```


### Galaxy in Docker container

```
docker run -d -p 8080:80 -p 8021:21 -p 8022:22 bgruening/galaxy-stable

docker commit <hash> bgruening/galaxy-stable-latest
```

docker run -v ~/WORK/MOBSuiteWrapper/galaxy_tools/tools:/export bgruening/planemo lint /export/filtering.xml

```
docker run  --privileged=true --rm -it -p 9001:9001 -p 9021:21 --name galaxy-stable bgruening/galaxy-stable

docker run -i -t -p 8080:80 bgruening/galaxy-stable /bin/bash
startup

docker run  -it -p 8080:80 -p 8022:22 --rm bgruening/galaxy-stable-latest

        
docker commit b97cd47c7535 bgruening/galaxy-stable-latest



CONTAINER ID        IMAGE                     COMMAND              CREATED             STATUS              PORTS                                                                     NAMES
41a56b18c866        bgruening/galaxy-stable   "/usr/bin/startup"   47 seconds ago      Up 51 seconds       443/tcp, 8800/tcp, 9002/tcp, 0.0.0.0:9021->21/tcp, 0.0.0.0:9080->80/tcp   galaxy-stable
k
```


supervisorctl stop  postgresql
supervisorctl stop  galaxy:*
supervisorctl stop nginx
supervisorctl stop cron


docker run --rm -it ubuntu:16.04






    sleep 5
        if [[ $NONUSE != *"postgres"* ]]
    then
        echo "Starting postgres"
        supervisorctl start postgresql
    fi
        wait_for_postgres
        if [[ $NONUSE != *"cron"* ]]
    then
        echo "Starting cron"
        supervisorctl start cron
    fi
            if [[ $NONUSE != *"proftp"* ]]
    then
        echo "Starting ProFTP"
        supervisorctl start proftpd
    fi
            if [[ $NONUSE != *"reports"* ]]
    then
        echo "Starting Galaxy reports webapp"
        supervisorctl start reports
    fi
            if [[ $NONUSE != *"nodejs"* ]]
    then
        echo "Starting nodejs"
        supervisorctl start galaxy:galaxy_nodejs_proxy
    fi
            if [[ $NONUSE != *"condor"* ]]
    then
        echo "Starting condor"
        supervisorctl start condor
             
        

###Build


https://github.com/dockerfile/ubuntu

```
docker build -f ./Dockerfile .  -t galaxy:latest

#run it
docker run -v /Users/kirill/WORK/MOBSuiteWrapper/galaxy_tools/tools/mob_suite/:/galaxy/tools/mobsuite -p 8080:8080 -it  2006f743dc5a  /bin/bash

#once inside the container
cd /galaxy && ./run.sh

#to stop the server
./run.sh stop

```

To install tools manually follow [this](https://galaxyproject.org/admin/tools/add-tool-tutorial/) instruction.

```
#publish docker image
docker login --username=kbessonov 
DOCKER_ID_USER="kbessonov"
docker tag galaxy $DOCKER_ID_USER/galaxy18.05

docker push $DOCKER_ID_USER/galaxy18.05
```

keys with issues 
--min_mob_evalue
--min_rpp_evalue
--min_rep_evalue
--min_con_evalue

```Traceback (most recent call last):
  File "/Users/kirill/miniconda/envs/mob_suite/bin/mob_recon", line 11, in <module>
    sys.exit(main())
  File "/Users/kirill/miniconda/envs/mob_suite/lib/python3.6/site-packages/mob_suite/mob_recon.py", line 410, in main
    if value > 1:
TypeError: '>' not supported between instances of 'str' and 'int'
```




mob_recon  --num_threads ${GALAXY_SLOTS:-4} --infile "MOB-Recon: Chromosomal sequences" --unicycler_contigs --run_circlator --run_typer --outdir '.'

import subprocess
input="/galaxy/database/files/000/dataset_25.dat"
subprocess.check_output( "head -n1 input | sed s/[\s+\|\:]//g" , shell=True)