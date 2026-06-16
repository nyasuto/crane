# src/crane/references_kuo.py
"""Kuo (2002) "Energetics of actively powered locomotion using the simplest
walking model", J. Biomech. Eng. 124(1):113-120 の公表知見。Phase 5a 検証ゲート。

URL: https://pubmed.ncbi.nlm.nih.gov/11871597/  取得日: 2026-06-17

対象モデル: Garcia/Kuo の simplest walking model（point-mass hip + 無質量脚 + 微小足質量）を
水平地面歩行へ拡張し、push-off（toe-off）インパルスと stance 脚トルクで駆動する powered
simplest walker。本 Phase 5a の実装（push-off インパルス P で COM を redirection し、
push-off 仕事 = P²/2）の直接の出典。

【原典本文は paywall（J. Biomech. Eng. / ASME）。本文 PDF は未取得。以下は (a) PubMed の
abstract を verbatim で確認した記述、(b) 明確に Kuo 2002 を引用する二次資料・索引本文から
確認した記述、を区別して記録する。Claude の記憶を数値ゲートにはしない（CLAUDE.md）。】

verbatim 確認（一次 abstract, PubMed 11871597, 2026-06-17 取得）:
    - "Powered walking was explored using an impulse applied at toe-off immediately
      before heel strike, and a torque applied on the stance leg."
    - "Although both methods can supply energy through mechanical work on the center
      of mass, the toe-off impulse is four times less costly because it decreases the
      collision loss at heel strike."  ← 有名な "1/4 the cost"（pre-emptive push-off が
      hip 仕事より 4 倍効率的）。これは abstract から直接確認できた。
    - "An idealized model yields a set of simple power laws relating the toe-off
      impulses and effective spring constant to the speed and step length of the
      corresponding gait."  ← push-off インパルスが speed/step length と power law で
      結ぶことは abstract に明記。ただし指数（2 乗等）は abstract には印刷されていない。

二次資料・索引本文から確認（指数の根拠。原典本文 PDF 未取得のため secondhand）:
    - Kuo, Donelan, Kram (2002) "Mechanical work for step-to-step transitions is a
      major determinant of the metabolic cost of human walking", J. Exp. Biol.
      205(23):3717 (https://journals.biologists.com/jeb/article/205/23/3717,
      2026-06-17 取得・abstract verbatim):
        "average negative and positive external mechanical work rates increased with
         the fourth power of step length (... r2=0.96)";
        "Metabolic rate also increased with the fourth power of step length (...
         r2=0.95), and linearly with mechanical work rate."  ← 固定 step frequency 下。
    - Kuo, Donelan, Ruina (2005) "Energetic consequences of walking like an inverted
      pendulum: step-to-step transitions", Exerc. Sport Sci. Rev. 33(2):88-97
      (https://pubmed.ncbi.nlm.nih.gov/15821430/, 2026-06-17): step-to-step transition
      の負仕事 W(-) ∝ |strike impulse|²、step-to-step 仕事は (walking speed × step
      length) の 2 乗で増加、と索引本文が記述。
    - 「正/負の外的力学仕事は push-off / collision インパルスの 2 乗に比例（Kuo 2002）」は
      上記群の索引本文・review が一貫して Kuo 2002 に帰属。点質量へのインパルス的仕事は
      解析的に ½·(impulse)²/m となり、本実装の push-off 仕事 = P²/2 と整合。

総括: push-off 仕事 ∝ push-off²（= インパルス²、本実装は厳密に P²/2）は abstract の power-law
記述＋二次資料で支持。力学コストが（固定 step length では）speed² で増える step-to-step
transition の関係も二次資料で支持。指数 "2"/"4" の verbatim 一次印刷は本文 paywall のため
未取得（abstract は指数を印刷しない）。よって以下のフラグは「実装と整合する定性関係が
読んだ資料で支持される」根拠で True とし、数値照合は None（捏造禁止）のままエネルギー収支
ゲートに委譲する。
"""

PROVENANCE: str = (
    "Kuo, A.D. (2002) 'Energetics of actively powered locomotion using the simplest "
    "walking model', J. Biomech. Eng. 124(1):113-120. 一次 abstract を "
    "https://pubmed.ncbi.nlm.nih.gov/11871597/ (取得日 2026-06-17) で verbatim 確認: "
    "'the toe-off impulse is four times less costly because it decreases the collision "
    "loss at heel strike'（pre-emptive push-off が hip 仕事より 4 倍効率的）, および "
    "'a set of simple power laws relating the toe-off impulses and effective spring "
    "constant to the speed and step length'. 原典本文 PDF は ASME paywall のため未取得。"
    "指数（push-off 仕事 ∝ push-off²、step-to-step 仕事 ∝ (speed×step length)²、仕事率 "
    "∝ step length⁴）は Kuo を明確に引用する二次資料で確認: Kuo,Donelan,Kram (2002) JEB "
    "205:3717 (https://journals.biologists.com/jeb/article/205/23/3717) と "
    "Kuo,Donelan,Ruina (2005) ESSR 33:88-97 (https://pubmed.ncbi.nlm.nih.gov/15821430/), "
    "いずれも 2026-06-17 取得。照合可能な単一数値は abstract に印刷されないため None とし、"
    "エネルギー収支（push-off 仕事 = P²/2）ゲートに委譲する（捏造禁止）。"
)

# Kuo の力学コスト ∝ 速度²、push-off 仕事 ∝ push-off²（本実装は厳密に P²/2）。
# True の根拠: push-off 仕事 ∝ push-off²（= インパルス²）は一次 abstract の power-law 記述
# ＋二次資料（W ∝ |impulse|²）で支持され、解析的にも点質量インパルス仕事 = ½P²/m と一致。
# 力学コストが固定 step length で speed² 増は二次資料（step-to-step transition）で支持。
# 指数の verbatim 一次印刷は本文 paywall で未取得である点は PROVENANCE に明記済み。
COT_SCALES_WITH_SPEED_SQUARED: bool = True

# 照合可能な数値（条件付き）。一次 abstract に単一の照合数値は印刷されていないため None。
# "4 倍効率（1/4 cost）" は定性比較であり単一の力学コスト数値ではないので PUBLISHED_VALUE には
# 載せず、PROVENANCE と下記 DESC に記録するにとどめる（数値ゲートには使わない）。
PUBLISHED_VALUE: float | None = None
PUBLISHED_VALUE_DESC: str = (
    "一次 abstract に照合用の単一数値印刷なし。定性結果: pre-emptive push-off は hip 仕事より "
    "4 倍効率的（toe-off impulse is four times less costly, Kuo 2002 abstract）。"
)
