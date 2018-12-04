# SigMA 

### Pre-Reqs: 
- [docker](https://docs.docker.com/engine/installation/)

### Running the app:
- `docker run -d -p 3242:3242 scottx611x/sigma`
- open http://localhost:3242
- Voila!

### Local Development:
- Make your edits to the application code
- Run `sudo sh build_image.sh` -> you should get a docker image ID in the console if things built properly
- Open http://localhost:3242 and poke around the app
- If you're satisfied with the changes create a pull request and/or notify @scottx611x

### Continuous integration/Deployment ToDos:
- [ ] New docker images will be built, tagged, and pushed to: [docker hub](https://hub.docker.com/r/scottx611x/sigma/) if a `PR` or `TAG` are noticed.
- We use [watchtower](https://github.com/v2tec/watchtower) to keep the current deployment up to date. Watchtower will detect new docker images that have been pushed, and update our remote docker gracefully.


### Deployment:
- AWS m5.large instance
- Accessible at: http://compbio.med.harvard.edu/sigMA
- `compbio.med.harvard.edu` providing http/ws proxy in `/etc/apache2/sites-available/compbio.med.harvard.edu.conf` to AWS instance
- Apache Proxy Config:
  ```
     Define sigMA_host <hostname>
     RewriteEngine on
     RewriteCond %{HTTP:Upgrade} =websocket
     RewriteRule /sigMA/(.*) ws://${sigMA_host}/$1 [P,L]
     RewriteCond %{HTTP:Upgrade} !=websocket
     RewriteRule /sigMA/(.*) http://${sigMA_host}/$1 [P,L]
     RewriteRule /sigMA http://compbio.med.harvard.edu/sigMA/
     ProxyPreserveHost On
     ProxyPass /sigMA http://${sigMA_host}/
     ProxyPassMatch ^/sigMA$ http://${sigMA_host}
     ProxyPassMatch ^/sigMA/(.*)$ http://${sigMA_host}/$1
     ProxyPassReverse /sigMA http://${sigMA_host}
     SetEnv force-proxy-request-1.0 1
     SetEnv proxy-nokeepalive 1
  ```  
- Docker container w/ run command: 
    + `sudo docker run --restart always -d -p 80:3242 -v /var/log/sigMA_logs/:/var/log/shiny-server/ scottx611x/sigma`
- Watchtower command:
    + `docker run -d --name watchtower -v /var/run/docker.sock:/var/run/docker.sock -e "REPO_USER=<dockerhub username>" -e "REPO_PASS=<dockerhub password" v2tec/watchtower --debug`
