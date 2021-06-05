---
title: "How to log in AWS"
---

## Today's IP

```sh
IP="54.188.62.211"
```

## Connect to pre-configured AWS machine

```sh
KEY="path/to/your/key.pem"
chmod 600 ${KEY}
#
USER= #e.g. user=user2
IP= # e.g. IP=255.255.255.255
#
ssh -i ${KEY} ${USER}@${IP}
```

E.g., for `user1`: 

```sh
KEY="c1.pem"
chmod 600 ${KEY}
USER=user1
IP=54.188.62.211
ssh -i ${KEY} ${USER}@${IP}
```
