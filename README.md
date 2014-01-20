## Centrifugal　Pump　Design　Tool
Centrifugal Pump design tool with EA/LOX engine
Only use for Sigle-stage certrifugal pump.
It is limited the use of approximate culcuration.
This script is New BSD License.
cf. Basic Design on Pumps by Hirohisa Takeda

## 超小型ロケット用遠心ポンプの設計ツール
単段の遠心ポンプの設計ツールです。超小型ロケットエンジン用を想定しているものです。
一般産業用の遠心ポンプの設計ではStepanoffの設計法が有名で一般的ですが、
ここでは電業社機械製作所の論文を参考に豊倉の設計法を採用しています。
修正BSDライセンスです。使用に関して一切の責任を負いかねます。

### 使用
pythonを使用しています。3次方程式の解を簡単にするためにscipyのfsolveを使用しました。
    $ python centrifugalpump.py
で出力されます。パラメータ変更ははソースファイルをいじって下さい。

### 決めるべき設計パラメータ
- ロケット推力、比推力、O/F
- ポンプ回転数
- 揚程H
- 流量Q
- 比速度Ns
- 羽根出口角beta2
- ポンプ効率
- 圧力低下係数λ
- 出口径D2
- 羽根枚数Z

### 参考文献
+ 武田祐久,"ポンプ設計の基礎",電業社機械,Vol.29,No.2,(2005),pp.7-11
http://www.dmw.co.jp/technical/pdf/no57.pdf
+ Stepanoff, A. J., “Centrifugal and Axial Flow Pumps(2nd Ed.),” John Wiley & Sons Inc., (1957).

+ 豊倉富太郎，武田裕久，“タ－ボポンプの新しい設計線図について”，タ－ボ機械第 32 回講演会，(1994), pp. 37-42.

+ 大嶋政夫,“遠心ポンプ羽根車に対するステパノフによる設計定数に関する一考察”，タ－ボ機械,Vol. 29, (2001), pp. 221-227

+ 渡辺啓悦,　吉田義樹,　"ロケットエンジン用ターボポンプ主羽根車小径化に関する検討",JAXAリポジトリ
http://repository.tksc.jaxa.jp/dr/prc/japan/contents/AA0063736000/63736000.pdf

+ 青木宏, 志村隆, 藁科彰吾, 上條謙二郎,"極低温上段エンジン用ターボポンプの設計および開発", 日本航空学会論文集,Vol.53, No.617,(2005), pp. 257-265
https://www.jstage.jst.go.jp/article/jjsass/53/617/53_617_257/_pdf