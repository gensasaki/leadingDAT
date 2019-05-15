library(KFAS)
library(ggplot2)
library(ggfortify)
library(zoo)

# 2003年1から2018年12月までのデータ
df <- read.csv("tourists_utf-8.csv", header=T, sep=",")
data <- ts(df$イタリア, start=c(2003, 1), frequency=12)
# log10_data <- ts(log10(df$中国), start=c(2003, 1), frequency=12)

ggplot(df$イタリア)
plot(data, type="l", ylab="The Number of Italian visited to Japan")
autoplot(data, xlab="Time", ylab="The Number of Italian visited to Japan")

# モデル定義 (ガウシアン)
# ローカルレベルモデル + 固定季節変動
mod1_g <- SSModel(data ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = 0), H = NA)
# ローカルレベルモデル + 可変季節変動
mod2_g <- SSModel(data ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = NA), H = NA)
# ローカル線形トレンドモデル + 固定季節変動
mod3_g <- SSModel(data ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = 0), H = NA)
# ローカル線形トレンドモデル + 可変季節変動
mod4_g <- SSModel(data ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = NA), H = NA)

# 未知パラメータの推定
fit1_g <- fitSSM(mod1_g, numeric(2), method = "BFGS")
fit2_g <- fitSSM(mod2_g, numeric(3), method = "BFGS")
fit3_g <- fitSSM(mod3_g, numeric(2), method = "BFGS")
fit4_g <- fitSSM(mod4_g, numeric(3), method = "BFGS")

# カルマンフィルタ・カルマンスムーザの実行
kfs1_g <- KFS(fit1_g$model)
kfs2_g <- KFS(fit2_g$model)
kfs3_g <- KFS(fit3_g$model)
kfs4_g <- KFS(fit4_g$model)
# 最大対数尤度
logLik1_g <- kfs1_g$logLik - sum(kfs1_g$Finf > 0) * log(2 * pi) / 2
logLik2_g <- kfs2_g$logLik - sum(kfs2_g$Finf > 0) * log(2 * pi) / 2
logLik3_g <- kfs3_g$logLik - sum(kfs3_g$Finf > 0) * log(2 * pi) / 2
logLik4_g <- kfs4_g$logLik - sum(kfs4_g$Finf > 0) * log(2 * pi) / 2
# AIC（赤池情報量規準）AIC = -2* 最大対数尤度 + 2*(未知パラメータ数 + 非定常な状態成分数)
AIC1_g <- -2 * logLik1_g + 2 * (2 + 12)
AIC2_g <- -2 * logLik2_g + 2 * (3 + 12)
AIC3_g <- -2 * logLik3_g + 2 * (2 + 13)
AIC4_g <- -2 * logLik4_g + 2 * (3 + 13)

# 原系列とトレンド（水準）成分の平滑化推定値および予測値
plot(data, type="l")
lines(kfs1_g$alphahat[, "level"], col=2)
lines(kfs2_g$alphahat[, "level"], col=3)
lines(kfs3_g$alphahat[, "level"], col=4)
lines(kfs4_g$alphahat[, "level"], col=5)






# モデル定義（ポアソン）
mod3_p <- SSModel(data ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = 0), distribution = "poisson")
mod4_p <- SSModel(data ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = NA), distribution = "poisson")

# 未知パラメータの推定
fit3_p <- fitSSM(mod3_p, c(-15, -10), nsim=1000, method = "BFGS")
fit4_p <- fitSSM(mod4_p, c(-15, -10), nsim=1000, method = "BFGS")

# カルマンフィルタ・カルマンスムーザの実行
kfs3_p <- KFS(fit3_p$model, nsim=1000)
kfs4_p <- KFS(fit4_p$model, nsim=1000)
# AIC
AIC3_p <- 2 * fit3_p$optim.out$value + 2 * (1 + 13)
AIC4_p <- 2 * fit4_p$optim.out$value + 2 * (2 + 13)

predict(fit3_p$model, interval="prediction", level=0.95, nsim=10000)
predict(fit4_p$model, interval="prediction", level=0.95, nsim=10000)



# 残差分析
mod0 <- SSModel(log(data) ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = NA), H=NA)
fit0 <- fitSSM(mod0, numeric(3))
kfs0 <- KFS(fit0$model, smoothing=c("state", "mean", "disturbance"))

autoplot(log(data))

autoplot(kfs0$model)

plot(log(data), type="l", ylab="The number of Italian visited to Japan (log)")
lines(kfs0$alphahat[, "level"], col=8)

plot(rstandard(kfs0, "pearson"), ylab="The number of Italian visited to Japan (log)")
abline(h=c(-1.96, 1.96), lty=3)

plot(rstandard(kfs0, "state")[, 1], ylab="The number of Italian visited to Japan (log)")
abline(h=c(-1.96, 1.96), lty=3)

# 異常な残差への対処
# 2011年3月、4月、5月、6月のデータを異常値として欠測値扱い
data[99] <- NA
data[100] <-NA
data[101] <- NA
data[102] <- NA
data[103] <- NA
data[104] <- NA

# 2016年3月以降の水準シフト干渉変数の定義
ShiftLevel <- c(rep(0, 158), rep(1, 22))

# 水準シフト干渉変数を加えたモデル
modShift <- SSModel(log(data) ~ SSMtrend(1, Q=NA)
                    + SSMseasonal(12, Q=NA, sea.type="dummy")
                    + ShiftLevel, H=NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing=c("state", "mean", "disturbance") )

plot(log(data), type="l", ylab="The number of Italian visited to Japan (log)")
lines(kfsShift$alphahat[, "level"], col=8)

plot(rstandard(kfsShift, "pearson"), ylab="The number of Italian visited to Japan (log)")
abline(h=c(-1.96, 1.96), lty=3)

plot(rstandard(kfsShift, "state")[, 1], at=2003:2017, ylab="The number of Italian visited to Japan (log)")
abline(h=c(-1.96, 1.96), lty=3)


# 2018年1〜10月の予測1
predict_data <- ts(c(data, rep(NA, 10)), start=c(2003, 1), frequency=12)

predict_mod2_g <- SSModel(predict_data ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = NA), H = NA)
predict_fit2_g <- fitSSM(predict_mod2_g, numeric(3), method = "BFGS")
predict_kfs2_g <- KFS(predict_fit2_g$model)
predict_logLik2_g <- predict_kfs2_g$logLik - sum(predict_kfs2_g$Finf > 0) * log(2 * pi) / 2
predict_AIC2_g <- -2 * predict_logLik2_g + 2 * (3 + 12)
preNA <- predict(predict_fit2_g$model, interval="prediction", level=0.95, nsim=10000)
confNA <- predict(predict_fit2_g$model, interval="confidence", level=0.95, nsim=10000)

plot(data, type="o", lty=3, xlab="Time", ylab="The Number of Italian visited to Japan")
lines(181:190, type="o", lty=3, col=8)

# 2018年1〜10月の予測2
# mod2_g <- SSModel(data ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = NA), H = NA)
# fit2_g <- fitSSM(mod2_g, numeric(3), method = "BFGS")
# kfs2_g <- KFS(fit2_g$model)
# logLik2_g <- kfs2_g$logLik - sum(kfs2_g$Finf > 0) * log(2 * pi) / 2
# AIC2_g <- -2 * logLik2_g + 2 * (3 + 12)
# preNA <- predict(predict_fit2_g$model, interval="prediction", level=0.95, nsim=10000, n.ahead=10)
# confNA <- predict(predict_fit2_g$model, interval="confidence", level=0.95, nsim=10000, n.ahead=10)

# plot(data, type="o", lty=3, xlab="Time", ylab="The Number of Italian visited to Japan")
# lines(181:190, data[181:190], type="o", lty=3, col=8)
# lines(181:190, confNA[, 1], lwd=2)
# lines(181:190, confNA[, 2])
# lines(181:190, confNA[, 3])
# lines(181:190, preNA[, 2], lty=2)
# lines(181:190, preNA[, 3], lty=2)


### カレンダー効果
# 各月の曜日集計
dates <- seq(as.Date("2003-01-01"), as.Date("2018-10-31"), by=1)
weeks <- table(substr(dates, 1, 7), weekdays(dates, T))
sun <- weeks[, "Sun"]
mon <- weeks[, "Mon"] - sun
tue <- weeks[, "Tue"] - sun
wed <- weeks[, "Wed"] - sun
thu <- weeks[, "Thu"] - sun
fry <- weeks[, "Fri"] - sun
sat <- weeks[, "Sat"] - sun
calendar <- cbind(mon, tue, wed, thu, fry, sat)

# うるう年2月のダミー変数
leapyear <- rownames(weeks) %in% c("2004-02", "2008-02", "2012-02", "2016-02")

# カレンダー効果（うるう年効果・曜日効果）のあるモデル
modCalendar <- SSModel(predict_data ~ SSMtrend(1, Q=NA) + SSMseasonal(12, Q=0) + leapyear + calendar, H=NA)
fitCalendar <- fitSSM(modCalendar, numeric(3), method="BFGS")
kfsCalendar <- KFS(fitCalendar$model)

# 対数尤度・AIC
logLikCalendar <- kfsCalendar$logLik - sum(kfsCalendar$Finf > 0) * log(2 * pi) / 2
AICCalendar <- -2 * logLikCalendar + 2*(2 + 20)
plot(kfsCalendar$muhat - kfsCalendar$alphahat[, "level"], type="l", ylab="The Number of Italian visited to Japan")
