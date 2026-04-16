# rwe_mama_datasus
 Estudo RWE — Trajetórias terapêuticas no SUS (Câncer de Mama, SP, 2019)
# RWE Câncer de Mama — Trajetórias Terapêuticas no SUS

Estudo de RWE usando dados públicos do DATASUS para comparar trajetórias terapêuticas no SUS-SP. 

---

## O que eu quis responder

A pergunta central é: Qual é a efetividade das trajetórias terapêuticas e como elas se refletem na 
persistência ao tratamento?

O framework PECO ficou assim:

- **P** — Mulheres ≥ 18 anos, CID C50.x, SP, 2019, new user design
- **E** — Quimioterapia + trastuzumabe (identifiquei pelo regex `HERCE|TRAST` nos campos de esquema terapêutico)
- **C** — Quimioterapia isolada, sem anti-HER2 registrado
- **O** — TTD: tempo até um gap > 90 dias entre APACs consecutivas

Por que TTD e não sobrevida? Tentei o linkage SIA × SIM. Não funcionou bem (detalhes abaixo). O TTD é um desfecho válido de qualquer forma — para Market Access, saber quanto tempo o paciente fica em tratamento no mundo real talvez seja até mais acionável do que sobrevida.

---

## Como rodar

R ≥ 4.3. Os notebooks precisam rodar em ordem:

```
00_setup.Rmd           → cria os diretórios
01_bronze.Rmd          → baixa SIA-AQ SP 2019 + geometria dos municípios
02_silver.Rmd          → limpeza, encoding, linkage SIM, variáveis geo
03_gold.Rmd            → monta a coorte analítica e os fluxos geográficos
04_analise_rwe.Rmd     → análise principal (descritiva, IPTW, KM, Cox)
4_a_aprofundamento.Rmd → modelo preditivo
4_b_aprofundamento.Rmd → geoprocessamento e Cox com interação
```

A ingestão usa `microdatasus::fetch_datasus()` — precisa de internet e pode demorar alguns minutos dependendo da conexão com o FTP do DATASUS.

```r
remotes::install_github("rfsaldanha/microdatasus")

install.packages(c("arrow", "sfarrow", "dplyr", "lubridate", "geobr", "sf",
                   "survival", "survminer", "WeightIt", "cobalt",
                   "gtsummary", "caret", "pROC", "ggstats"))
```

---

## O pipeline

Organizei em três camadas. Nada muito sofisticado, mas garante que eu sempre consiga voltar ao dado bruto se precisar.

**Bronze** — dado bruto, sem tocar. SIA-AQ SP 2019, SIM-DO SP 2019 e geometria municipal do geobr. Tudo em `.parquet`.

Dois missings que importam logo de cara:
- `AQ_ESTADI` (estadiamento): 12,7% ausente
- `AQ_ESQU_P2`: 21,7% ausente — o campo de esquema terapêutico é texto livre legado, às vezes o nome do medicamento cabe no P1, às vezes transborda pro P2. Isso provavelmente subestima o grupo trastuzumabe.

**Silver** — limpeza. Os clássicos do DATASUS: encoding latin1 → UTF-8 (o `iconv` é obrigatório nesses campos de texto livre), missings disfarçados como `"99"`, `"000000"`, datas como string que precisam virar `Date`. Criei a variável de deslocamento intermunicipal (`AP_MUNPCN ≠ AP_UFMUN`) e calculei os centroides municipais via geobr para o mapa.

Sobre o linkage SIA × SIM: tentei linkage probabilístico por município de residência + faixa etária (±1 ano) + CID C50 + janela de 30 dias após a última APAC. A taxa não convergiu para nada plausível epidemiologicamente. Faz sentido: os dados são desidentificados, São Paulo tem municípios enormes, e a combinação de atributos disponíveis não é específica o suficiente para distinguir pacientes. Desisti do linkage e pivotei para TTD como desfecho principal. Os óbitos flagrados dentro do próprio SIA (`AP_OBITO`) continuam sendo usados.

**Gold** — coorte analítica. New user design (só quem iniciou tratamento em 2019), exposição definida por Intention-to-Treat (presença de trastuzumabe nas 3 primeiras competências para mitigar immortal time bias), time zero na primeira `dt_inicio` real — não na competência de faturamento, que tem ruído de calendário administrativo.

---

## O que encontrei

Antes dos números, um aviso importante sobre interpretação — voltarei a isso nas limitações.

**IPTW:** usei `WeightIt` com regressão logística (ATE). Apareceram pesos extremos (máximo 73,5), o que indica que alguns pacientes tinham probabilidade muito baixa de pertencer ao grupo oposto — provavelmente reflexo do confundimento por subtipo molecular, que não está nos dados. Fiz trimming no percentil 99,5. O Love Plot mostrou balanceamento adequado (SMD < 0,10) após o ajuste.

**Kaplan-Meier ponderado por IPTW** — as curvas se separam cedo e a diferença é estatisticamente significativa.

**Cox duplamente robusto** (pesos IPTW + `robust = TRUE`):

| | HR | IC 95% | p |
|---|---|---|---|
| Trastuzumabe vs. Quimio isolada | 1,48 | 1,16 – 1,88 | 0,002 |
| CACON vs. UNACON | 0,65 | — | 0,002 |
| Deslocamento intermunicipal | 0,77 | — | 0,067 |

O teste de Schoenfeld mostrou violação marginal da proporcionalidade em `deslocamento` (p = 0,013) e `idade`. A variável de exposição principal cumpriu o pressuposto, então a estimativa do HR do trastuzumabe não está comprometida. Optei por não introduzir termos de interação com o tempo para essas variáveis — adicionaria complexidade sem ganho na estimativa que interessa.

---

## Aprofundamento A — modelo preditivo

Regressão logística com validação cruzada 10-fold para identificar pacientes com alto risco de gap > 90 dias. A ideia prática é: esse modelo poderia alimentar a priorização de um PSP — quem acionar primeiro, com que recurso.

Features: tratamento, idade, estadiamento, tipo de serviço, interação `tipo_serv × deslocamento`. Avaliação por AUC-ROC honesta (out-of-fold).

---

## Aprofundamento B — geoprocessamento

O HR de deslocamento com tendência protetora (0,77) levantou uma hipótese: o efeito provavelmente não é da distância em si, mas de para onde o paciente está indo. Testei isso com um modelo Cox com interação `tipo_serv * deslocamento`:

| Estrato | HR |
|---|---|
| CACON, mesmo município | 0,41 (p < 0,001) |
| UNACON, deslocamento intermunicipal | 0,51 (p < 0,001) |
| **CACON + deslocamento intermunicipal** | **2,71 (p < 0,001)** |

Tratar em CACON perto de casa é muito protetor. Viajar para uma UNACON melhor também é (seletividade — quem viaja é o paciente motivado). Mas viajar para um CACON anula tudo isso e mais um pouco. A logística da alta complexidade é o gargalo real.

O mapa de desire lines materializa isso — linhas de origem (município de residência) → destino (CACON), espessura proporcional ao volume de APACs. Limitação: trabalho no centroide municipal, não no endereço do estabelecimento. Para Market Access regional é suficiente.

---

## Limitações — a parte que mais importa

**A limitação mais grave da análise não é o linkage SIM. É a ambiguidade do próprio desfecho.**

Um gap > 90 dias pode ser:
- Abandono por toxicidade ou dificuldade de acesso — o que queremos capturar
- **Alta clínica** — paciente completou o protocolo e foi dispensada
- Óbito não capturado
- Mudança de serviço

No contexto do câncer de mama isso é especialmente problemático. O trastuzumabe adjuvante tem **duração protocolar definida de 12 meses**. Uma fração relevante do grupo exposto vai ter um gap > 90 dias simplesmente porque o tratamento terminou conforme planejado. O SIA não tem campo confiável para "alta por conclusão de protocolo" — o `AP_MOTSAI` existe mas é mal preenchido.

Isso cria um viés direcional: o grupo trastuzumabe acumula "eventos" de TTD artificiais porque o protocolo HER2+ tem fim programado, enquanto alguns esquemas de quimioterapia isolada continuam indefinidamente no contexto paliativo. O HR > 1 encontrado é consistente com esse viés — e provavelmente está sendo parcialmente explicado por ele, não só pelo subtipo molecular.

O número de 1,48 não significa que trastuzumabe leva ao abandono mais rápido. Significa que trajetórias HER2+ terminam mais cedo no SIA — por qualquer razão. A interpretação correta é "trajetórias terapêuticas distintas no SUS", não "eficácia diferente".

Outras limitações relevantes:

- **Subtipo molecular não observado** — HER2+ e HER2- são doenças biologicamente distintas. O IPTW balanceia o que está nos dados, mas o principal fator de confundimento não está nos dados.
- **AQ_ESQU_P2 com 21,7% missing** — subestima o grupo exposto. Viés conservador, mas existe.
- **Schoenfeld violado em deslocamento** — HR do deslocamento instável no tempo; não compromete a exposição principal.
- **Interação CACON + deslocamento pode capturar gravidade** — pacientes que viajam para CACONs podem ser casos mais graves, não apenas mais distantes. Confundimento residual não controlável com esses dados.
- **Recorte 2019** — deliberado para evitar viés COVID, mas limita generalização.

---

## Decisões técnicas rápidas


IPTW em vez de PSM porque não descarta observações. Com amostras que já têm limitações de covariáveis, perder linhas por falta de par no matching parece pior do que lidar com pesos extremos via trimming.

