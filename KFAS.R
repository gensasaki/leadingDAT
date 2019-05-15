### KFAS

library(KFAS)

setwd("~/Work/R_projects/leadingDAT_LB1/")

Weight <- ts(read.delim("data/Weight.dat", sep = " ", header = FALSE))
Bodyfat <- ts(read.delim ("data/Bodyfat.dat", sep = " ", header = FALSE))

# ローカルレベルモデルの定義
mod <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)
# 未知パラメータQ, Hの推定（最尤法）
fit <- fitSSM(mod, numeric(2), method = "BFGS")
# カルマンフィルタ・カルマンスムーザの実行
kfs <- KFS(fit$model)


## 平滑化と欠測値の補完
# 平滑化状態の信頼区間
alphahatconf <- predict(fit$model, interval = "confidence", level = 0.95)
# 欠測値の補完
WeightNA <- Weight[c(1:20, rep(NA, 20), 41:60)] # 21~40日目をNAに置き換え
modNA <- SSModel(WeightNA ~ SSMtrend(1, Q = NA), H = NA)
fitNA <- fitSSM(modNA, c(0, 0), method = "BFGS")
confNA <- predict(fitNA$model, interval = "confidence", level = 0.95)
preNA <- predict(fitNA$model, interval = "prediction", level = 0.95)


## 長期予測
# 長期予測
mod50 <- SSModel(Weight[1:50] ~ SSMtrend(1, Q = NA), H = NA)
fit50 <- fitSSM(mod50, numeric(2), method = "BFGS")
conf50 <- predict(fit50$model, interval = "confidence", n.ahead = 10)
pre50 <- predict(fit50$model, interval = "prediction", n.ahead = 10)

# 長期予測（予測期間を欠測値NAとして予測するやり方）
Weight50 <- Weight[c(1:50, rep(NA, 10))] # 51日目以降をNA（欠損）に置き換え
mod50NA <- SSModel(Weight50 ~ SSMtrend(1, Q = NA), H = NA)
fit50NA <- fitSSM(mod50NA, numeric(2), method = "BFGS")
conf50NA <- predict(fit50NA$model, interval = "confidence", level = 0.95)
pre50NA <- predict(fit50NA$model, interval = "prediction", level = 0.95)


## 2変量ローカルレベルモデルの実装と推定結果
modSUTSE <- SSModel(cbind(Weight, Bodyfat) ~ SSMtrend(1, Q = matrix(NA, 2, 2)), H = matrix(NA, 2, 2))
fitSUTSE <- fitSSM(modSUTSE, numeric(6), method = "BFGS")
kfsSUTSE <- KFS(fitSUTSE$model)


## 平滑化トレンドモデルの実装
# 水準・傾きを持つモデル
# 1.ローカルレベルモデル
# 2.平滑化トレンドモデル
# 3.ローカル線型トレンドモデル

# モデル定義
mod1 <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)
mod2 <- SSModel(Weight ~ SSMtrend(2, Q = list(0, NA)), H = NA)
mod3 <- SSModel(Weight ~ SSMtrend(2, Q = list(NA, NA)), H = NA)

# 未知パラメータの推定
fit1 <- fitSSM(mod1, numeric(2), method = "BFGS")
fit2 <- fitSSM(mod2, numeric(2), method = "BFGS")
fit3 <- fitSSM(mod3, numeric(3), method = "BFGS")

# カルマンフィルタ・カルマンスムーザの実行
kfs1 <- KFS(fit1$model)
kfs2 <- KFS(fit2$model)
kfs3 <- KFS(fit3$model)


## 対数尤度とAICによるモデル選択
# 対数尤度
logLik1 <- kfs1$logLik - sum(kfs1$Finf > 0) * log(2 * pi) / 2
logLik2 <- kfs2$logLik - sum(kfs2$Finf > 0) * log(2 * pi) / 2
logLik3 <- kfs3$logLik - sum(kfs3$Finf > 0) * log(2 * pi) / 2

# AIC
AIC1 <- -2 * logLik1 + 2 * (2 + 1)
AIC2 <- -2 * logLik1 + 2 * (2 + 2)
AIC3 <- -2 * logLik1 + 2 * (3 + 2)


## 基本構造時系列モデルの実装
sales <- read.csv("data/sales.csv", fileEncoding="Shift_JIS") # データ（csv形式）の読み込み
fabric <- ts(sales[, 2]) # 織物衣服業の月次販売額データ
fabricNA <- c(fabric[1:120], rep(NA, 24)) # 直近2年間をNA（欠測=予測対象）
# モデル定義（平滑化トレンド(mod3, mod4)は省略）
mod1 <- SSModel(fabric ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = 0), H = NA)
mod2 <- SSModel(fabricNA ~ SSMtrend(1, Q = NA) + SSMseasonal(12, Q = NA), H = NA)
mod3 <- SSModel(fabricNA ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = 0), H = NA)
mod4 <- SSModel(fabricNA ~ SSMtrend(2, Q = list(0, NA)) + SSMseasonal(12, Q = NA), H = NA)
# 未知パラメータの推定
fit1 <- fitSSM(mod1, numeric(2), method = "BFGS")
fit2 <- fitSSM(mod2, numeric(3), method = "BFGS")
fit3 <- fitSSM(mod3, numeric(2), method = "BFGS")
fit4 <- fitSSM(mod4, numeric(3), method = "BFGS")
# カルマンフィルタ・カルマンスムーザの実行
kfs1 <- KFS(fit1$model)
kfs2 <- KFS(fit2$model)
kfs3 <- KFS(fit3$model)
kfs4 <- KFS(fit4$model)
# 最大対数尤度
logLik1 <- kfs1$logLik - sum(kfs1$Finf > 0) * log(2 * pi) / 2
logLik2 <- kfs2$logLik - sum(kfs2$Finf > 0) * log(2 * pi) / 2
logLik3 <- kfs3$logLik - sum(kfs3$Finf > 0) * log(2 * pi) / 2
logLik4 <- kfs4$logLik - sum(kfs4$Finf > 0) * log(2 * pi) / 2
# AIC（赤池情報量規準）
AIC1 <- -2 * logLik1 + 2 * (2 + 12)
AIC2 <- -2 * logLik2 + 2 * (3 + 12)
AIC3 <- -2 * logLik3 + 2 * (2 + 13)
AIC4 <- -2 * logLik4 + 2 * (3 + 13)
# 図の描画（横軸と縦線の描画, 軸の設定など一部省略）
# (a)原系列とトレンド（水準）成分の平滑化推定値および予測値
plot(fabric, type = "l")
lines(kfs2$alphahat[, "level"], col = 3)
lines(kfs4$alphahat[, "level"], col = 4)
# (b)季節成分の平滑化推定値および予測値
plot(kfs2$alphahat[, "sea_dummy1"], type = "l", col = 3)
lines(kfs4$alphahat[, "sea_dummy1"], col = 4)
# (c)平滑化後の残差と予測誤差
plot(fabric - kfs2$muhat, type = "l", col = 3)
lines(fabric - kfs4$muhat, col = 4)


## カレンダー効果の実装
# 各月の曜日集計
dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = 1)
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
leapyear <- rownames(weeks) %in% c("2004--2", "2008-02", "2012-02")
# カレンダー効果（うるう年効果・曜日効果）のあるモデル
modCalendar <- SSModel(fabricNA ~ SSMtrend(2, Q = list(0, NA))
                       + SSMseasonal(12, sea.type = "dummy", Q = 0) + leapyear + calendar, H = NA)
fitCalendar <- fitSSM(modCalendar, numeric(2), method = "BFGS")
kfsCalendar <- KFS(fitCalendar$model)


## カレンダーの解析例
# 対数尤度・AIC
logLikCalendar <- kfsCalendar$logLik - sum(kfsCalendar$Finf > 0) * log(2 * pi) / 2
AICCalendar <- -2 * logLikCalendar + 2 * (2 + 10)
plot(kfsCalendar$muhat - kfsCalendar$alphahat[, "level"], type = "l")


## 状態空間モデルの残差分析の例
machine <- ts(sales$Machinery...Equipment)
mod0 <- SSModel(log(machine) ~ SSMtrend(1, Q = NA)
                + SSMseasonal(12, Q = NA, sea.type = "dummy"), H = NA)
fit0 <- fitSSM(mod0, numeric(3))
kfs0 <- KFS(fit0$model, smoothing = c("state", "mean", "disturbance"))
# 図の描画
plot(log(machine), type = "l")
lines(kfs0$alphahat[, "level"], col = 8)
plot(rstandard(kfs0, "pearson", ylim = c(-6, 6)))
abline(h = c(-1.96, 1.96), lty = 3)
plot(rstandard(kfs0, "state")[, 1], ylim = c(-6, 6))
abline(h = c(-1.96, 1.96), lty = 3)


## 異常な残差への対処
# 2010年11月, 12月のデータを異常値として除外（NAとして欠測値扱い）
machineNA <- machine
machineNA[sales[, 1] %in% c("2010年11月", "2010年12月")] <- NA
# 2011年8月以降の水準シフト干渉変数の定義
ShiftLevel <- (1:nrow(sales) >= which(sales[, 1] == "2011年8月"))
# 水準シフト干渉変数を加えたモデル
modShift <- SSModel(log(machineNA) ~ SSMtrend(1, Q = NA)
                    + SSMseasonal(12, Q = NA, sea.type = "dummy") + ShiftLevel, H = NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing = c("state", "mean", "disturbance"))


## 月次死者数の解析：モデルの実装
# 状態空間モデルの定義
fire <- ts(read.delim("data/fire.dat", sep = " ", header = FALSE))
modPois <- SSModel(fire ~ SSMtrend(2, Q = list(0, NA))
                   + SSMseasonal(12, Q = NA), distribution = "poisson")
# 最尤法による未知パラメータの推定
fitPois <- fitSSM(modPois, c(-15, -10), nsim = 1000, method = "BFGS")
# 状態推定
kfsPois <- KFS(fitPois$model, nsim = 1000)
# 長期予測（予測区間）
prePoisNA <- predict(fitPois$model, interval = "prediction", level = 0.95, nsim = 10000)

